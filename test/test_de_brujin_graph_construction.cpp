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
void testDeBruijnGraph(const mxx::comm& comm, const std::string & filename, const std::string test) {

	NodeMapType idx(comm);

	TIMER_INIT(test);

	if (comm.rank() == 0)
		INFOF("RANK %d / %d: Testing %s", comm.rank(), comm.size(), test.c_str());

	TIMER_START(test);
	idx.template build<SeqParser>(filename, comm);
	TIMER_END(test, "build", idx.local_size());

	TIMER_START(test);
	auto query = readForQuery<NodeMapType>(filename, comm);
	TIMER_END(test, "read query", query.size());

	// for testing, query 1% (else could run out of memory.  if a kmer exists r times, then we may need r^2/p total storage.
	if (idx.local_size() > 1000) {
		TIMER_START(test);
		unsigned seed = comm.rank() * 23;
		sample(query, query.size() / 100, seed);
		TIMER_END(test, "select 1%", query.size());
	}

	// process query
	// query
	TIMER_START(test);
	auto results = idx.find(query);
	TIMER_END(test, "query", results.size());

  for (auto result : results) {
    std::cout << result << std::endl;
  }

	TIMER_REPORT_MPI(test, comm.rank(), comm);
}

using Alphabet = bliss::common::DNA;
using KmerType = bliss::common::Kmer<21, Alphabet, WordType>;
using EdgeEncoder = bliss::common::DNA16;

template <typename EdgeEnc>
using CountNodeMapType = bliss::de_bruijn::de_bruijn_nodes_distributed<
		KmerType, bliss::de_bruijn::node::edge_counts<EdgeEnc, int32_t>, int,
		bliss::kmer::transform::lex_less,
		bliss::kmer::hash::farm>;

template <typename EdgeEnc>
using ExistNodeMapType = bliss::de_bruijn::de_bruijn_nodes_distributed<
    KmerType, bliss::de_bruijn::node::edge_exists<EdgeEnc>, int,
    bliss::kmer::transform::lex_less,
    bliss::kmer::hash::farm>;
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

	std::string filename("/home/tpan/src/bliss/test/data/test.debruijn.small.fastq");
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



	::std::cerr<<"Using DNA16 to present each edge" << ::std::endl;
	testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<CountNodeMapType>, bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, count."));

//
//	::std::cerr<<"Using ASCII to present each edge" << ::std::endl;
//	testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine_ascii<CountNodeMapType>, bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, count."));


  ::std::cerr<<"Using DNA16 to represent each edge" << ::std::endl;
  testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine<ExistNodeMapType>,  bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, existence."));

//
//  ::std::cerr<<"Using ASCII to represent each edge" << ::std::endl;
//  testDeBruijnGraph< bliss::de_bruijn::de_bruijn_engine_ascii<ExistNodeMapType >,  bliss::io::FASTQParser >(comm, filename, ::std::string("ST, hash, dbg construction, existence."));


	/*synchronization*/
	MPI_Barrier(comm);

	//////////////  clean up MPI.
	MPI_Finalize();

	//INFOF("M Rank %d called MPI_Finalize", rank);

	return 0;
}

