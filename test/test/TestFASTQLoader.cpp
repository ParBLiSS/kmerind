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
      FILE *fp = fopen(fileName.c_str(), "r");
      fseek(fp, offset , SEEK_SET);
      size_t read = fread_unlocked(result, 1, length, fp);
      fclose(fp);

      if (read < 0) throw std::logic_error("ERROR: fread_unlocked failed.");
    }




  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <template <typename> class SeqParser>
  size_t read_file(const std::string & filename, const mxx::comm & _comm) {


//    constexpr size_t overlap = KP::kmer_type::size;  specify as 0 - which allows overlap to be computed.

    // TODO: check if prefetch makes a difference.
    using FileLoaderType = bliss::io::FileLoader<CharType, 0, SeqParser, false, true>; // raw data type :  use CharType

    //====  now process the file, one L1 block (block partition by MPI Rank) at a time

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
      unsigned char * gold = new unsigned char[len + 1];

      readFilePOSIX(filename, offset, len, gold);

      bool same = ::std::equal(partition.begin(), partition.end(), gold);

      if (!same) std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (fileloader)! first chars " << *(partition.begin()) << " gold " << (*gold)
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
  template <template <typename> class SeqParser>
  static size_t read_file_mpi_subcomm(const std::string & filename, const mxx::comm& _comm) {

    // split the communcator so 1 proc from each host does the read, then redistribute.
    mxx::comm group = _comm.split_shared();
    mxx::comm group_leaders = _comm.split(group.rank() == 0);

//    constexpr size_t overlap = KP::kmer_type::size;  specify as 0 - which allows overlap to be computed.

    // raw data type :  use CharType.   block partition at L1 and L2.  no buffering at all, since we will be copying data to the group members anyways.
    using FileLoaderType = bliss::io::FileLoader<CharType, 0, SeqParser, false, true,
        bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >,  bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >>;
    using L2BlockType = bliss::io::DataBlock<unsigned char*, ::bliss::partition::range<size_t>, bliss::io::NoBuffer>;

    typename FileLoaderType::RangeType file_range;
    L2BlockType block;

    TIMER_INIT(file);
    {

      typename FileLoaderType::RangeType range;
      ::std::vector<CharType> data;

      ::std::vector<typename FileLoaderType::RangeType> ranges;
      ::std::vector<size_t> send_counts;
      typename FileLoaderType::L1BlockType partition;

      // first load the file using the group loader's communicator.
      TIMER_START(file);
      if (group.rank() == 0) {  // ensure file loader is closed properly.
        //==== create file Loader. this handle is alive through the entire building process.
        FileLoaderType loader(filename, group_leaders, group.size());  // for member of group_leaders, each create g_size L2blocks.

        // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
        partition = loader.getNextL1Block();

        //====  now compute the send counts and ranges to be scattered.
        ranges.resize(group.size());
        send_counts.resize(group.size());
        for (int i = 0; i < group.size(); ++i) {
          ranges[i] = loader.getNextL2Block(i).getRange();
          send_counts[i] = ranges[i].size();
        }

        // scatter the data .  this call here relies on FileLoader still having the memory mapped.
        // TODO; use iterators!
        data = mxx::scatterv(&(*partition.begin()), send_counts, 0, group);

        // send the file range.
        file_range = loader.getFileRange();

      } else {
        // replicated here because the root's copy needs to be inside the if clause - requires FileLoader to be open for the sender.

        // scatter the data.  this is the receiver end of the thing.
        data = mxx::scatterv(&(*partition.begin()), send_counts, 0, group);

      }

      // send the file range to rest of group
      mxx::datatype range_dt = mxx::get_datatype<typename FileLoaderType::RangeType >();
      MPI_Bcast(&file_range, 1, range_dt.type(), 0, group);


      // scatter the ranges to rest of group
      // TODO: mxx::bcast function
      range = mxx::scatter_one(ranges, 0, group);
      // now create the L2Blocks from the data  (reuse block)
      block.assign(&(data[0]), &(data[0]) + range.size(), range);
      TIMER_END(file, "open", data.size());


      // check partition is same
      size_t len = block.getRange().size();
      size_t offset = block.getRange().start;
      unsigned char * gold = new unsigned char[len + 1];

      readFilePOSIX(filename, offset, len, gold);

      bool same = ::std::equal(block.begin(), block.end(), gold);

      if (!same) std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (subcomm)! first chars " << *(partition.begin()) << " gold " << (*gold)
              << " range " << partition.getRange() << std::endl;

      delete [] gold;



//      // not reusing the SeqParser in loader.  instead, reinitializing one.
//      TIMER_START(file);
//      SeqParser<typename L2BlockType::iterator> l2parser;
//      l2parser.init_parser(block.begin(), file_range, range, range, comm);
//      TIMER_END(file, "mark_seqs", est_size);

    }

    TIMER_REPORT_MPI(file, _comm.rank(), _comm);

    return block.getRange().size();
  }



  /**
   * @brief read a file's content and generate kmers, place in a vector as return result.
   * @note  static so can be used wihtout instantiating a internal map.
   * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.   template template parameter, param is iterator
   * @tparam KmerParser   parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
   */
  template <template <typename> class SeqParser, size_t overlap = 0>
  static size_t read_file_direct(const std::string & filename, const mxx::comm & _comm) {


    //====  now process the file, one L1 block (block partition by MPI Rank) at a time

    // TODO: check if prefetch makes a difference.  specify overlap as 0 - which allows overlap to be computed.
    ::bliss::io::MemData partition;
    partition.valid_range.start = 0;
    partition.valid_range.end = 0;


    TIMER_INIT(file);
    {  // ensure that fileloader is closed at the end.

      TIMER_START(file);
      //==== create file Loader
      partition = ::bliss::io::parallel::memmap::load_file<SeqParser<unsigned char*> >(filename, _comm, overlap);
      TIMER_END(file, "open", partition.getRange().size());

      // check partition is same
      size_t len = partition.getRange().size();
      size_t offset = partition.getRange().start;
      unsigned char * gold = new unsigned char[len + 1];

      readFilePOSIX(filename, offset, len, gold);

      bool same = ::std::equal(partition.begin(), partition.end(), gold);

      if (!same) std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (direct)! first chars " << *(partition.begin()) << " gold " << (*gold)
          << " range " << partition.getRange() << std::endl;

      delete [] gold;


//      // not reusing the SeqParser in loader.  instead, reinitializing one.
//      TIMER_START(file);
//      SeqParser<bliss::io::DataType*> l1parser;
//      l1parser.init_parser(partition.data, partition.parent_range, partition.mem_range, partition.valid_range, _comm);
//      TIMER_END(file, "mark_seqs", est_size);

    }

    TIMER_REPORT_MPI(file, _comm.rank(), _comm);
    return partition.valid_range.size();
  }



template <template <typename> class SeqParser>
void testIndex(const mxx::comm& comm, const std::string & filename, std::string test ) {

  TIMER_INIT(test);

  if (comm.rank() == 0) BL_INFOF("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  TIMER_START(test);
  size_t result = read_file<SeqParser>(filename, comm);
  TIMER_END(test, "build", result);




  TIMER_REPORT_MPI(test, comm.rank(), comm);

}


template <template <typename> class SeqParser>
void testIndexSubComm(const mxx::comm& comm, const std::string & filename, std::string test ) {

  TIMER_INIT(test);

  if (comm.rank() == 0) BL_INFOF("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  TIMER_START(test);
  size_t result = read_file_mpi_subcomm<SeqParser>(filename, comm);
  TIMER_END(test, "build", result);




  TIMER_REPORT_MPI(test, comm.rank(), comm);

}


template <template <typename> class SeqParser, typename KmerType>
void testIndexDirect(const mxx::comm& comm, const std::string & filename, std::string test ) {

  TIMER_INIT(test);

  if (comm.rank() == 0) BL_INFOF("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

  TIMER_START(test);
  size_t result = read_file_direct<SeqParser, KmerType::size>(filename, comm);
  TIMER_END(test, "build", result);




  TIMER_REPORT_MPI(test, comm.rank(), comm);

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
      filename.append("/test/data/test.fastq");

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
  
  if (which != -1) std::cin.get();


  using Alphabet = bliss::common::DNA;
  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;



  if (which == -1 || which == 1)
  {
  testIndex<bliss::io::FASTQParser > (comm, filename, "ST, hash, all read, count index.");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 2)
  {
  testIndexSubComm< bliss::io::FASTQParser > (comm, filename, "ST, hash, read_and_dist, count index.");
    MPI_Barrier(comm);
}

  if (which == -1 || which == 3)
  {
  testIndexDirect<bliss::io::FASTQParser , KmerType> (comm, filename, "ST, hash, read_and_dist, count index.");
    MPI_Barrier(comm);
}

  //////////////  clean up MPI.
  MPI_Finalize();

  //BL_INFOF("M Rank %d called MPI_Finalize", rank);



  return 0;
}
