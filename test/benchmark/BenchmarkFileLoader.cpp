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
 * @file    test_threads.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *

 */

#include "bliss-config.hpp"

#include <unistd.h>  // get hostname


#include <functional>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <chrono>
#include <iostream>  // for system("pause");
#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"
#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "index/quality_score_iterator.hpp"

#include "index/kmer_index.hpp"

#include "utils/timer.hpp"



void readFilePOSIX(const std::string &fileName, const size_t offset,
						  const size_t length, unsigned char* result)
{
  FILE *fp = fopen64(fileName.c_str(), "r");
  fseeko64(fp, offset , SEEK_SET);
  size_t read = fread_unlocked(result, 1, length, fp);
  fclose(fp);

  if (read < 0) throw std::logic_error("ERROR: fread_unlocked failed.");
}


#if defined(USE_FASTQ_PARSER)
#define PARSER_TYPE ::bliss::io::FASTQParser
#elif defined(USE_FASTA_PARSER)
#define PARSER_TYPE ::bliss::io::FASTAParser
#endif



  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <template <typename> class SeqParser, bool prefetch = false>
  size_t read_file(const std::string & filename, const mxx::comm & _comm) {


//    constexpr size_t overlap = KP::kmer_type::size;  specify as 0 - which allows overlap to be computed.

    // TODO: check if prefetch makes a difference.
    using FileLoaderType = bliss::io::FileLoader<CharType, 0, SeqParser, false, prefetch>; // raw data type :  use CharType

    //====  now process the file, one L1 block (block partition by MPI Rank) at a time
    if (_comm.rank() == 0) printf("using prefetch? %s\n", (prefetch ? "y" : "n"));

    typename FileLoaderType::L1BlockType partition;
    TIMER_INIT(file);
    {  // ensure that fileloader is closed at the end.

      TIMER_START(file);
      //==== create file Loader
      FileLoaderType loader(filename, _comm, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
      partition = loader.getNextL1Block();
      TIMER_END(file, "open", partition.getRange().size());


      // check partition is same
      size_t len = partition.getRange().size();
      size_t offset = partition.getRange().start;

      TIMER_START(file);
      unsigned char * gold = new unsigned char[len + 1];
      readFilePOSIX(filename, offset, len, gold);
      TIMER_END(file, "posix", len);


      TIMER_START(file);
      auto diff_iters = ::std::mismatch(partition.begin(), partition.end(), gold);
      TIMER_END(file, "compare", len);


      if (static_cast<size_t>(std::distance(partition.begin(), partition.end())) != len)
        std::cout << "ERROR: block size is not same as range size." << std::endl;

      if (diff_iters.first != partition.end())
        std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (subcomm)! diff at offset " << (offset + std::distance(partition.begin(), diff_iters.first))
            << " val " << *(diff_iters.first) << " gold " << *(diff_iters.second)
            << " range " << partition.getRange() << std::endl;


      delete [] gold;

//      // not reusing the SeqParser in loader.  instead, reinitializing one.
//      TIMER_START(file);
//      auto l1parser = loader.getSeqParser();
//      l1parser.init_parser(partition.begin(), loader.getFileRange(), partition.getRange(), partition.getRange(), _comm);
//      TIMER_END(file, "mark_seqs", est_size);

    }


    TIMER_REPORT_MPI(file, _comm.rank(), _comm);
    return partition.getRange().size();
  }


  /**
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.  template template parameter, param is iterator
   */
  template <template <typename> class SeqParser, bool prefetch = false>
  static size_t read_file_mpi_subcomm(const std::string & filename, const mxx::comm& _comm) {

    // split the communcator so 1 proc from each host does the read, then redistribute.
    mxx::comm group = _comm.split_shared();
    mxx::comm group_leaders = _comm.split(group.rank() == 0);

//    constexpr size_t overlap = KP::kmer_type::size;  specify as 0 - which allows overlap to be computed.
    if (_comm.rank() == 0) printf("using prefetch? %s\n", (prefetch ? "y" : "n"));

    // raw data type :  use CharType.   block partition at L1 and L2.  no buffering at all, since we will be copying data to the group members anyways.
    using FileLoaderType = bliss::io::FileLoader<CharType, 0, SeqParser, false, prefetch,
        bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >,  bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >>;
    using L2BlockType = bliss::io::DataBlock<unsigned char*, ::bliss::partition::range<size_t>, bliss::io::NoBuffer>;

    typename FileLoaderType::RangeType file_range;
    L2BlockType block;

    TIMER_INIT(file_subcomm);
    {

      typename FileLoaderType::RangeType range;
      ::std::vector<CharType> data;

      ::std::vector<typename FileLoaderType::RangeType> ranges;
      ::std::vector<size_t> send_counts;
      typename FileLoaderType::L1BlockType partition;

      // first load the file using the group loader's communicator.
      if (group.rank() == 0) {  // ensure file loader is closed properly.
          TIMER_START(file_subcomm);
        //==== create file Loader. this handle is alive through the entire building process.
        FileLoaderType loader(filename, group_leaders, group.size());  // for member of group_leaders, each create g_size L2blocks.

        // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
        partition = loader.getNextL1Block();
        TIMER_END(file_subcomm, "L1 block", data.size());

        //====  now compute the send counts and ranges to be scattered.
        TIMER_START(file_subcomm);
        ranges.resize(group.size());
        send_counts.resize(group.size());
        for (int i = 0; i < group.size(); ++i) {
          auto b = loader.getNextL2Block(i);
          ranges[i] = b.getRange();
          send_counts[i] = ranges[i].size();
        }
        TIMER_END(file_subcomm, "L2 range", data.size());

        TIMER_START(file_subcomm);
        // scatter the data .  this call here relies on FileLoader still having the memory mapped.
        data = mxx::scatterv(&(*partition.begin()), send_counts, 0, group);
        TIMER_END(file_subcomm, "scatter", data.size());

        // send the file range.
        file_range = loader.getFileRange();

      } else {
        // replicated here because the root's copy needs to be inside the if clause - requires FileLoader to be open for the sender.


          TIMER_START(file_subcomm);
          TIMER_END(file_subcomm, "L1 block", data.size());
          TIMER_START(file_subcomm);
          TIMER_END(file_subcomm, "L2 range", data.size());

        // scatter the data.  this is the receiver end of the thing.
          TIMER_START(file_subcomm);
        data = mxx::scatterv(&(*partition.begin()), send_counts, 0, group);
        TIMER_END(file_subcomm, "scatter", data.size());


      }
      if (data[0] != '@') std::cout << "rank " << _comm.rank() << " mxx data begins with " << data[0] << std::endl;

      TIMER_START(file_subcomm);

      // send the file range to rest of group
      mxx::datatype range_dt = mxx::get_datatype<typename FileLoaderType::RangeType >();
      MPI_Bcast(&file_range, 1, range_dt.type(), 0, group);


      // scatter the ranges to rest of group
      // TODO: mxx::bcast function
      range = mxx::scatter_one(ranges, 0, group);
      // now create the L2Blocks from the data  (reuse block)
      block.assign(&(data[0]), &(data[0]) + range.size(), range);
      TIMER_END(file_subcomm, "createL2", data.size());


      // check partition is same
      size_t len = block.getRange().size();
      size_t offset = block.getRange().start;

      TIMER_START(file_subcomm);
      unsigned char * gold = new unsigned char[len + 1];
      readFilePOSIX(filename, offset, len, gold);
      TIMER_END(file_subcomm, "posix", len);


      TIMER_START(file_subcomm);
      auto diff_iters = ::std::mismatch(block.begin(), block.end(), gold);
      TIMER_END(file_subcomm, "compare", len);
//
//      if (std::distance(block.begin(), block.end()) != len)
//        std::cout << "ERROR: block size is not same as range size." << std::endl;

      if (diff_iters.first != block.end())
        std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (subcomm)! diff at offset " << (offset + std::distance(block.begin(), diff_iters.first))
            << " val " << *(diff_iters.first) << " gold " << *(diff_iters.second)
            << " range " << block.getRange() << std::endl;

      delete [] gold;



//      // not reusing the SeqParser in loader.  instead, reinitializing one.
//      TIMER_START(file);
//      SeqParser<typename L2BlockType::iterator> l2parser;
//      l2parser.init_parser(block.begin(), file_range, range, range, comm);
//      TIMER_END(file, "mark_seqs", est_size);

    }

    TIMER_REPORT_MPI(file_subcomm, _comm.rank(), _comm);

    return block.getRange().size();
  }



  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <typename FileLoader, size_t overlap = 0>
  static size_t read_file_mpi_direct(const std::string & filename, const mxx::comm & _comm) {


    //====  now process the file, one L1 block (block partition by MPI Rank) at a time

    ::bliss::io::file_data partition;
    partition.valid_range_bytes.start = 0;
    partition.valid_range_bytes.end = 0;


    TIMER_INIT(file_direct);
    {  // ensure that fileloader is closed at the end.

      TIMER_START(file_direct);
      //==== create file Loader
      FileLoader loader(filename, overlap, _comm);
      TIMER_END(file_direct, "open", partition.getRange().size());

      TIMER_START(file_direct);
      //==== create file Loader
      partition = loader.read_file();
      TIMER_END(file_direct, "load", partition.getRange().size());


      // check partition is same
      size_t len = partition.getRange().size();
      size_t offset = partition.getRange().start;

      TIMER_START(file_direct);
      unsigned char * gold = new unsigned char[len + 1];
      readFilePOSIX(filename, offset, len, gold);
      TIMER_END(file_direct, "posix", len);


      TIMER_START(file_direct);
      auto diff_iters = ::std::mismatch(partition.begin(), partition.end(), gold);
      TIMER_END(file_direct, "compare", len);


      if (partition.data.size() != partition.in_mem_range_bytes.size())
        std::cout << "ERROR rank " << _comm.rank() << ": block size " << partition.data.size() << " is not same as range size "
        << partition.valid_range_bytes.size()
        << " in mem range size is " << partition.in_mem_range_bytes.size() << std::endl;

      if (diff_iters.first != partition.end())
        std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (subcomm)! start at " << offset 
		<< " diff at offset " << (offset + std::distance(partition.begin(), diff_iters.first))
            << " val " << *(diff_iters.first) << " gold " << *(diff_iters.second)
            << " range " << partition.getRange() << std::endl;


      delete [] gold;


//      // not reusing the SeqParser in loader.  instead, reinitializing one.
//      TIMER_START(file);
//      SeqParser<bliss::io::DataType*> l1parser;
//      l1parser.init_parser(partition.data, partition.parent_range, partition.mem_range, partition.valid_range, _comm);
//      TIMER_END(file, "mark_seqs", est_size);

    }

    TIMER_REPORT_MPI(file_direct, _comm.rank(), _comm);
    return partition.valid_range_bytes.size();
  }



template <template <typename> class SeqParser, bool prefetch>
void testIndex(const mxx::comm& comm, const std::string & filename, std::string test ) {


  if (comm.rank() == 0) printf("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  read_file<SeqParser, prefetch>(filename, comm);
}


template <template <typename> class SeqParser, bool prefetch>
void testIndexSubComm(const mxx::comm& comm, const std::string & filename, std::string test ) {


  if (comm.rank() == 0) printf("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  read_file_mpi_subcomm<SeqParser, prefetch>(filename, comm);
}


template <typename FileLoader, typename KmerType>
void testIndexDirect(const mxx::comm& comm, const std::string & filename, std::string test ) {

  if (comm.rank() == 0) printf("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  read_file_mpi_direct<FileLoader, KmerType::size>(filename, comm);
}

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {

  //////////////// init logging
  LOG_INIT();

  //////////////// parse parameters

  std::string filename;
      filename.assign(PROJ_SRC_DIR);

#if defined(USE_FASTQ_PARSER)
      filename.append("/test/data/test.fastq");
#elif defined(USE_FASTA_PARSER)
      filename.append("/test/data/test.fasta");
#endif

  if (argc > 1)
  {
    filename.assign(argv[1]);
  }

  int which = -1;
  if (argc > 2)
	which = atoi(argv[2]);


  int rank = 0;
	int nthreads = 1;
  //////////////// initialize MPI and openMP
#ifdef USE_MPI

  if (nthreads > 1) {

    int provided;

    // one thread will be making all MPI calls.
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    if (provided < MPI_THREAD_FUNNELED) {
      BL_ERRORF("The MPI Library Does not have thread support.");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);

  MPI_Barrier(comm);

  if (rank == 0) BL_INFOF("USE_MPI is set");
#else
  static_assert(false, "MPI used although compilation is not set to use MPI");
#endif
  
  //if (which != -1) std::cin.get();


  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;



  if (which == -1 || which == 1)
  {
  testIndex<PARSER_TYPE, true > (comm, filename, "file loader.");
    MPI_Barrier(comm);
  }

  if (which == -1 || which == 2)
  {
  testIndex<PARSER_TYPE, false > (comm, filename, "file loader.");
    MPI_Barrier(comm);
  }


  if (which == -1 || which == 3)
  {
  testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::mmap_file, PARSER_TYPE<unsigned char *> >, KmerType> (comm, filename, "mpi mmap");
    MPI_Barrier(comm);
  }

  if (which == -1 || which == 4)
  {
	  testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::stdio_file, PARSER_TYPE<unsigned char *> >, KmerType> (comm, filename, "mpi stdio");
    MPI_Barrier(comm);
  }
  if (which == -1 || which == 5)
  {
	  testIndexDirect<bliss::io::parallel::mpiio_file<PARSER_TYPE<unsigned char *> >, KmerType> (comm, filename, "mpi-io");
    MPI_Barrier(comm);
  }


  if (which == -1 || which == 6)
  {
  testIndexSubComm< PARSER_TYPE, true > (comm, filename, "fileloader with mpi subcomm.");
    MPI_Barrier(comm);
  }

  if (which == -1 || which == 7)
  {
  testIndexSubComm< PARSER_TYPE, false > (comm, filename, "fileloader with mpi subcomm.");
    MPI_Barrier(comm);
  }

  //////////////  clean up MPI.
  MPI_Finalize();

  //BL_INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
