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

#include <string>
#include <sstream>
#include <iostream>  // for system("pause");

#include "utils/logging.h"

#include "io/file.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/exception_handling.hpp"

#include "tclap/CmdLine.h"

#include "utils/mxx_fast_comm.hpp"


template <typename ITER>
bool validate(const std::string &fileName, const size_t offset,
						  const size_t length, ITER const & first, ITER const & last)
{
  FILE *fp = fopen64(fileName.c_str(), "r");
  fseeko64(fp, offset, SEEK_SET);

  if (static_cast<size_t>(std::distance(first, last)) != length)
    std::cout << "ERROR: block size is not same as range size." << std::endl;


  using valtype = typename std::iterator_traits<ITER>::value_type;
  constexpr size_t tmp_size = 1024 * 1024;

  valtype tmp[tmp_size];
  size_t count = 0;
  ITER iter = first;
  size_t read_size = tmp_size;

  for (size_t l = 0; l < length; l += tmp_size) {
    read_size = std::min(tmp_size, length - l);

	  count = fread_unlocked(tmp, sizeof(valtype), read_size, fp);  // pointer is advanced.


		if (count < read_size) {
		  if (feof(fp)) {
		    // no problem here
		  }  else if (ferror(fp)) {
		    printf("ERROR during read.  read %ld, less than requested size %ld\n", count, read_size);
		  }
		} else if (count > read_size) {
		  printf("ERROR during read.  read %ld, more than requested size %ld\n", count, read_size);
		}
	  auto diff_iters = ::std::mismatch(tmp, tmp + count, iter);
	  iter += count;

	  if (diff_iters.first != (tmp + count)) {
		  std::cout << "ERROR: diff at offset " << (offset + l + std::distance(tmp, diff_iters.first))
		  	  << " internal offset " << std::distance(tmp, diff_iters.first)
		  	  << " val " << *(diff_iters.second) << " gold " << *(diff_iters.first) << " count " << count << std::endl;
		  return false;
	  }
  }

  fclose(fp);

  return (iter == last);
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
  template <typename FileLoader, size_t overlap = 0>
  static size_t read_file_mpi_direct(const std::string & filename, const mxx::comm & _comm) {


    //====  now process the file, one L1 block (block partition by MPI Rank) at a time

    ::bliss::io::file_data partition;
    partition.valid_range_bytes.start = 0;
    partition.valid_range_bytes.end = 0;



    BL_BENCH_INIT(file_direct);
    {  // ensure that fileloader is closed at the end.


      BL_BENCH_START(file_direct);
      //==== create file Loader
      FileLoader loader(filename, overlap, _comm);
      BL_BENCH_END(file_direct, "open", partition.getRange().size());


      BL_BENCH_START(file_direct);
      //==== create file Loader
      partition = loader.read_file();
      BL_BENCH_END(file_direct, "load", partition.getRange().size());



      // check partition is same
      size_t len = partition.getRange().size();
      size_t offset = partition.getRange().start;


      if (partition.data.size() != partition.in_mem_range_bytes.size())
        std::cout << "ERROR rank " << _comm.rank() << ": block size " << partition.data.size() << " is not same as range size "
        << partition.valid_range_bytes.size()
        << " in mem range size is " << partition.in_mem_range_bytes.size() << std::endl;


      BL_BENCH_START(file_direct);
      bool same = validate(filename, offset, len, partition.begin(), partition.end());
      BL_BENCH_END(file_direct, "compare", len);

      if (!same) {
    	  std::cout << "ERROR: rank " << _comm.rank() << " NOT SAME (subcomm)! "
          	  << " range " << partition.getRange() << std::endl;
    	  exit(1);
      }

//      // not reusing the SeqParser in loader.  instead, reinitializing one.
//      BL_BENCH_START(file);
//      SeqParser<bliss::io::DataType*> l1parser;
//      l1parser.init_parser(partition.data, partition.parent_range, partition.mem_range, partition.valid_range, _comm);
//      BL_BENCH_END(file, "mark_seqs", est_size);

    }

    BL_BENCH_REPORT_MPI(file_direct, _comm.rank(), _comm);


    return partition.valid_range_bytes.size();
  }



//template <template <typename> class SeqParser, bool prefetch>
//void testIndex(const mxx::comm& comm, const std::string & filename, std::string test ) {
//
//
//  if (comm.rank() == 0) printf("RANK %d / %d: Testing %s\n", comm.rank(), comm.size(), test.c_str());
//
//  read_file<SeqParser, prefetch>(filename, comm);
//}
//
//
//template <template <typename> class SeqParser, bool prefetch>
//void testIndexSubComm(const mxx::comm& comm, const std::string & filename, std::string test ) {
//
//
//  if (comm.rank() == 0) printf("RANK %d / %d: Testing %s\n", comm.rank(), comm.size(), test.c_str());
//
//  read_file_mpi_subcomm<SeqParser, prefetch>(filename, comm);
//}
//

template <typename FileLoader, typename KmerType>
void testIndexDirect(const mxx::comm& comm, const std::string & filename, std::string test ) {

  if (comm.rank() == 0) printf("RANK %d / %d: Testing %s\n", comm.rank(), comm.size(), test.c_str());

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

  //////////////// initialize MPI and openMP

  mxx::env e(argc, argv);
  mxx::comm comm_world;
  comm_world.barrier();

  //////////////// parse parameters

  std::string filename;
  filename.assign(PROJ_SRC_DIR);
#if defined(USE_FASTQ_PARSER)
      filename.append("/test/data/test.fastq");
#elif defined(USE_FASTA_PARSER)
      filename.append("/test/data/test.fasta");
#endif

  int which = -1;
  int nnodes = -1;


  // Wrap everything in a try block.  Do this every time,
  // because exceptions will be thrown for problems.
  try {

    // Define the command line object, and insert a message
    // that describes the program. The "Command description message"
    // is printed last in the help text. The second argument is the
    // delimiter (usually space) and the last one is the version number.
    // The CmdLine object parses the argv array based on the Arg objects
    // that it contains.
    TCLAP::CmdLine cmd("Benchmark parallel file loading", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-n Bishop".
    TCLAP::ValueArg<std::string> fileArg("F", "file", "FASTQ file path", false, filename, "string", cmd);

    TCLAP::ValueArg<int> algoArg("A", "algo", "Algorithm id.  If absent, all.", false, -1, "int", cmd);
    TCLAP::ValueArg<int> nnodesArg("N", "nnodes", "Number of target nodes.  If absent, all.", false, comm_world.size(), "int", cmd);


    // Parse the argv array.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    filename = fileArg.getValue();
    which = algoArg.getValue();
    nnodes = nnodesArg.getValue();

    // Do what you intend.

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(-1);
  }

  // if nnodes is too large, then set to max.
  nnodes = std::min(comm_world.size(), nnodes);

  ::mxx::comm comm = comm_world.copy();

  bool is_fast = true;
  if (nnodes < comm_world.size()) {
    std::cout << "NOTE: comm_world is being split by bandwidth." << std::endl;
	  is_fast = ::bliss::mxx::get_fast_nodes(comm_world, e, nnodes);
	  comm = comm_world.split(is_fast ? 1 : 0);
  }

  // only run on fast nodes.
  if (is_fast) {
  
	  using Alphabet = bliss::common::DNA;
	  using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;


	  if (which == -1 || which == 5)
	  {
		  testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::mmap_file, PARSER_TYPE<typename ::bliss::io::file_data::const_iterator> >, KmerType> (comm, filename, "mpi mmap");
		  comm.barrier();
	  }

	  if (which == -1 || which == 6)
	  {
		  testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::stdio_file, PARSER_TYPE<typename ::bliss::io::file_data::const_iterator> >, KmerType> (comm, filename, "mpi stdio");
		  comm.barrier();
	  }

    if (which == -1 || which == 7)
    {
      testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::posix_file, PARSER_TYPE<typename ::bliss::io::file_data::const_iterator> >, KmerType> (comm, filename, "mpi posix");
      comm.barrier();
    }

    if (which == -1 || which == 8)
    {
      testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::mmap_file, PARSER_TYPE<typename ::bliss::io::file_data::const_iterator>, bliss::io::parallel::base_shared_fd_file >, KmerType> (comm, filename, "mpi fd mmap");
      comm.barrier();
    }

    if (which == -1 || which == 9)
    {
      testIndexDirect<bliss::io::parallel::partitioned_file<bliss::io::posix_file, PARSER_TYPE<typename ::bliss::io::file_data::const_iterator>, bliss::io::parallel::base_shared_fd_file >, KmerType> (comm, filename, "mpi fd posix");
      comm.barrier();
    }


	  if (which == -1 || which == 10)
	  {
		  testIndexDirect<bliss::io::parallel::mpiio_file<PARSER_TYPE<typename ::bliss::io::file_data::const_iterator> >, KmerType> (comm, filename, "mpi-io");
		  comm.barrier();
	  }



  }

  //BL_INFOF("M Rank %d called MPI_Finalize", rank);
  comm_world.barrier();

  //////////////  clean up MPI.
  //MPI_Finalize();


  return 0;
}
