/*
 * test_de_brujin_graph_construction.cpp
 *
 *  Created on: Aug 17, 2015
 *      Author: yongchao
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

#include <iterators/edge_iterator.hpp>
#include "wip/kmer_index.hpp"
#include "wip/de_bruijn_construct_engine.hpp"
#include "wip/de_bruijn_nodes_distributed.hpp"

#include "utils/timer.hpp"



template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  std::string extension = ::bliss::utils::file::get_file_extension(filename);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  if (extension.compare("fastq") == 0) {
    // default to including quality score iterators.
    IndexType::template read_file<::bliss::io::FASTQParser, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);
  } else {
    throw std::invalid_argument("input filename extension is not supported.");
  }

  return query;
}


template <typename IndexType, typename KmerType = typename IndexType::KmerType>
std::vector<KmerType> readForQuery_subcomm(const std::string & filename, MPI_Comm comm) {

  ::std::vector<KmerType> query;

  std::string extension = ::bliss::utils::file::get_file_extension(filename);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  if (extension.compare("fastq") == 0) {
    // default to including quality score iterators.
    IndexType::template read_file_mpi_subcomm<::bliss::io::FASTQParser, ::bliss::index::kmer::KmerParser<KmerType> >(filename, query, comm);
  } else {
    throw std::invalid_argument("input filename extension is not supported.");
  }

  return query;
}



template<typename KmerType>
void sample(std::vector<KmerType> &query, size_t n, unsigned int seed) {
	std::shuffle(query.begin(), query.begin() + ::std::min(4 * n, query.size()),
			std::default_random_engine(seed));
	query.erase(query.begin() + n, query.end());
}

template<typename NodeMapType, template <typename> class SeqParser>
void testDeBruijnGraph(MPI_Comm comm, const std::string & filename, const std::string test) {

	int nprocs = 1;
	int rank = 0;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &rank);

	NodeMapType idx(comm, nprocs);

	TIMER_INIT(test);

	if (rank == 0)
		INFOF("RANK %d / %d: Testing %s", rank, nprocs, test.c_str());

	TIMER_START(test);
	idx.template build<SeqParser>(filename, comm);
	TIMER_END(test, "build", idx.local_size());

	TIMER_START(test);
	auto query = readForQuery<NodeMapType>(filename, comm);
	TIMER_END(test, "read query", query.size());

	// for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
	TIMER_START(test);
	unsigned seed = rank * 23;
	sample(query, query.size() / 100, seed);
	TIMER_END(test, "select 1%", query.size());

	auto query_orig = query;

	auto query1 = query_orig;
	query1.resize(1);

	// query 1
	TIMER_START(test);
	auto results3 = idx.find(query1);
	TIMER_END(test, "query 1", results3.size());

	for (auto result : results3) {
	  std::cout << result << std::endl;
	}

	query1 = query_orig;
	query1.resize(1);

	query = query_orig;

	// process query
	// query
	TIMER_START(test);
	auto results = idx.find(query);
	TIMER_END(test, "query 1%", results.size());

  for (auto result : results) {
    std::cout << result << std::endl;
  }


	query = query_orig;

	TIMER_REPORT_MPI(test, rank, comm);
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

	std::string filename("/home/tpan/src/bliss/test/data/test.medium.fastq");
	if (argc > 1) {
		filename.assign(argv[1]);
	}

	cerr << "filename: " << filename << endl;

	int rank = 0;
	int size = 0;
	//////////////// initialize MPI and openMP
#ifdef USE_MPI

	/*initialize MPI environment*/
	MPI_Init(&argc, &argv);

	/*set the communicator; all processes participate in the processing*/
	MPI_Comm comm = MPI_COMM_WORLD;

	/*get the rank*/
	MPI_Comm_rank(comm, &rank);

	{
		char hostname[256];
		memset(hostname, 0, 256);
		gethostname(hostname, 256);
		//INFOF("Rank %d hostname [%s]", rank, hostname);
	}
	MPI_Comm_size(comm, &size);
	MPI_Barrier(comm);

	if (rank == 0)
		INFOF("USE_MPI is set");
#else
	static_assert(false, "MPI used although compilation is not set to use MPI");
#endif

	using Alphabet = bliss::common::DNA;
	using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;

	using Uint16NodeMapType = bliss::de_bruijn::de_bruijn_nodes_distributed<
			KmerType, bliss::de_bruijn::node::edge_counts<Alphabet, int32_t>, int,
			bliss::kmer::transform::lex_less,
			bliss::kmer::hash::farm>;

	using Uint8NodeMapType = bliss::de_bruijn::de_bruijn_nodes_distributed<
			KmerType, bliss::de_bruijn::node::edge_counts<Alphabet, int32_t>, int,
			bliss::kmer::transform::lex_less,
			bliss::kmer::hash::farm>;

	::std::cerr<<"Using uint16_t to present each edge" << ::std::endl;
	testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<Uint16NodeMapType>, bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, count."));

	::std::cerr<<"Using uint8_t to represent each edge" << ::std::endl;
	testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<Uint8NodeMapType>,  bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, count."));


  using Uint16NodeMapType2 = bliss::de_bruijn::de_bruijn_nodes_distributed<
      KmerType, bliss::de_bruijn::node::edge_exists<Alphabet>, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm>;

  using Uint8NodeMapType2 = bliss::de_bruijn::de_bruijn_nodes_distributed<
      KmerType, bliss::de_bruijn::node::edge_exists<Alphabet>, int,
      bliss::kmer::transform::lex_less,
      bliss::kmer::hash::farm>;

  ::std::cerr<<"Using 1bit to present each edge" << ::std::endl;
  testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<Uint16NodeMapType2>, bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, existence."));

  ::std::cerr<<"Using 1bit to represent each edge" << ::std::endl;
  testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<Uint8NodeMapType2>,  bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, existence."));

	/*synchronization*/
	MPI_Barrier(comm);

	//////////////  clean up MPI.
	MPI_Finalize();

	//INFOF("M Rank %d called MPI_Finalize", rank);

	return 0;
}

