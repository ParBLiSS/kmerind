/**
 * @file    kmer_index.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief	High level kmer indexing API
 * @details	3 primary Kmer Index classes are currently provided:
 * 			Kmer Count Index,
 * 			Kmer Position Index, and
 * 			Kmer Position + Quality score index.
 *
 *
 *		Currently, KmerIndex classes only supports FASTQ files, but the intent is to
 *		  allow replaceable module for file reading.
 *
 *		Data distribution is via replaceable Hash functions.  currently,
 *		  PrefixHasher is used to distribute to MPI processes based on the first log_2(P) bits of the kmer
 *		  InfixHasher is used to distribute to threads within a process, based on the next log_2(p) bits of the kmer
 *		  SuffixHasher is used for threadlocal unordered_map hashing
 *
 *		All are deterministic to allow simple lookup processes.
 *
 *
 * Copyright (c) 2014 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef KMERINDEX2_HPP_
#define KMERINDEX2_HPP_

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

#include "io/fastq_loader.hpp"
//#include "io/fasta_loader.hpp"
//#include "io/fasta_iterator.hpp"

#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"
#include <containers/distributed_unordered_map.hpp>
#include <containers/distributed_sorted_map.hpp>
// way too slow.  also not updated.  #include <containers/distributed_map.hpp>

#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "iterators/constant_iterator.hpp"
#include "index/quality_score_iterator.hpp"
#include "wip/kmer_hash.hpp"

#include "utils/timer.hpp"
#include "utils/string_utils.hpp"

namespace bliss
{
namespace index
{
namespace kmer
{

/**
 * @tparam MapType  	container type
 * @tparam KmerParser		functor to generate kmer (tuple) from input.  specified here so we specialize for different index.  note KmerParser needs to be supplied with a data type.
 */
template <typename MapType, typename KmerParser>
class Index {
protected:

	MapType map;

	MPI_Comm comm;
	int commSize;
	int commRank;

public:
	using KmerType = typename MapType::key_type;
	using IdType   = typename MapType::mapped_type;
	using TupleType =      std::pair<KmerType, IdType>;
	using Alphabet = typename KmerType::KmerAlphabet;

	Index(MPI_Comm _comm, int _comm_size) : map(_comm, _comm_size), comm(_comm) {
		MPI_Comm_size(_comm, &commSize);
		MPI_Comm_rank(_comm, &commRank);
	}

	virtual ~Index() {};

	MapType & get_map() {
		return map;
	}

	/**
	 * @brief  generate kmers or kmer tuples for 1 block of raw data.
	 *
	 * @tparam SeqParser		parser type for extracting sequences.  supports FASTQ and FASTA.
	 * @tparam BlockType		input partition type, supports in memory (vector) vs memmapped.
	 * @param partition
	 * @param result        output vector.  should be pre allocated.
	 */
	template <typename SeqParserType, typename KP, typename BlockType>
	static size_t read_block(BlockType & partition, std::vector<typename KP::value_type>& result) {

		// from FileLoader type, get the block iter type and range type
		using BlockIterType = typename BlockType::iterator;

		using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, SeqParserType>;
		using SeqType = typename ::std::iterator_traits<SeqIterType>::value_type;

		//== sequence parser type
		KP kmer_parser;

		//== process the chunk of data
		SeqType read;

		//==  and wrap the chunk inside an iterator that emits Reads.
		SeqIterType seqs_start(SeqParserType(), partition.begin(), partition.end(), partition.getRange().start);
		SeqIterType seqs_end(partition.end());

    ::fsc::back_emplace_iterator<std::vector<typename KP::value_type> > emplace_iter(result);

    size_t before = result.size();

		//== loop over the reads
		for (; seqs_start != seqs_end; ++seqs_start)
		{
			emplace_iter = kmer_parser(*seqs_start, emplace_iter);
		}

		return result.size() - before;
	}



	/**
	 * @tparam KmerParser		parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
	 * @tparam SeqParser		parser type for extracting sequences.  supports FASTQ and FASTA.
	 * @tparam T	data element type after reading from file.  May not be the same as Map's data type.
	 */
	template <typename SeqParser, typename KP = KmerParser>
	static size_t read_file(const std::string & filename, std::vector<typename KP::value_type>& result, MPI_Comm _comm) {

		int p, rank;
		MPI_Comm_size(_comm, &p);
		MPI_Comm_rank(_comm, &rank);

		using FileLoaderType = bliss::io::FileLoader<CharType, SeqParser, false, false>; // raw data type :  use CharType

		//====  now process the file, one L1 block (block partition by MPI Rank) at a time

		size_t before = result.size();

		TIMER_INIT(file);
		{  // ensure that fileloader is closed at the end.

			TIMER_START(file);
			//==== create file Loader
			FileLoaderType loader(_comm, filename, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
			typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
			TIMER_END(file, "open", partition.getRange().size());

			//== reserve
			TIMER_START(file);
			// modifying the local index directly here causes a thread safety issue, since callback thread is already running.
			// index reserve internally sends a message to itself.

			// call after getting first L1Block to ensure that file is loaded.
			size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + p - 1) / p;
			result.reserve(est_size);
			TIMER_END(file, "reserve", est_size);


			TIMER_START(file);
			//=== copy into array
			while (partition.getRange().size() > 0) {

				read_block<SeqParser, KP>(partition, result);

				partition = loader.getNextL1Block();
			}
			TIMER_END(file, "read", result.size());
		}

		TIMER_REPORT_MPI(file, rank, _comm);
		return result.size() - before;
	}


	/**
	 * @tparam KmerParser		parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
	 * @tparam SeqParser		parser type for extracting sequences.  supports FASTQ and FASTA.
	 * @tparam T	data element type after reading from file.  May not be the same as Map's data type.
	 */
	template <typename SeqParser, typename KP = KmerParser>
	static size_t read_file_mpi_subcomm(const std::string & filename, std::vector<typename KP::value_type>& result, MPI_Comm _comm) {

		int p, rank;
		MPI_Comm_size(_comm, &p);
		MPI_Comm_rank(_comm, &rank);

		// split the communcator so 1 proc from each host does the read, then redistribute.
		MPI_Comm group_leaders;
		MPI_Comm group;
		::std::tie(group_leaders, group) = ::mxx2::split_communicator_by_host(_comm);

		int g_size, g_rank;
		MPI_Comm_rank(group, &g_rank);
		MPI_Comm_size(group, &g_size);

		int gl_size = MPI_UNDEFINED;
		int gl_rank = MPI_UNDEFINED;
		if (group_leaders != MPI_COMM_NULL) {
			MPI_Comm_rank(group_leaders, &gl_rank);
			MPI_Comm_size(group_leaders, &gl_size);
		}

		// raw data type :  use CharType.   block partition at L1 and L2.  no buffering at all, since we will be copying data to the group members anyways.
		using FileLoaderType = bliss::io::FileLoader<CharType, SeqParser, false, false, bliss::partition::BlockPartitioner<bliss::partition::range<size_t> > >;

		size_t before = result.size();

		TIMER_INIT(file);
		{
			typename FileLoaderType::L2BlockType block;
			typename FileLoaderType::RangeType range;
			::std::vector<CharType> data;

			::std::vector<typename FileLoaderType::RangeType> ranges;
			::std::vector<size_t> send_counts;
			typename FileLoaderType::L1BlockType partition;

			// first load the file using the group loader's communicator.
			TIMER_START(file);
			size_t est_size = 1;
			if (group_leaders != MPI_COMM_NULL) {  // ensure file loader is closed properly.
				//==== create file Loader. this handle is alive through the entire building process.
				FileLoaderType loader(group_leaders, filename, g_size);  // for member of group_leaders, each create g_size L2blocks.

				// modifying the local index directly here causes a thread safety issue, since callback thread is already running.

				partition = loader.getNextL1Block();

				// call after getting first L1Block to ensure that file is loaded.
				est_size = (loader.getKmerCountEstimate(KmerType::size) + p - 1) / p;

				//====  now compute the send counts and ranges to be scattered.
				ranges.resize(g_size);
				send_counts.resize(g_size);
				for (size_t i = 0; i < g_size; ++i) {
					block = loader.getNextL2Block(i);

					ranges[i] = block.getRange();
					send_counts[i] = ranges[i].size();
				}

				// scatter the data .  this call here relies on FileLoader still having the memory mapped.
				data = ::mxx2::scatterv(partition.begin(), partition.end(), send_counts, group, 0);

			} else {
				// replicated here because the root's copy needs to be inside the if clause - requires FileLoader to be open for the sender.

				// scatter the data .  this call here relies on FileLoader still having the memory mapped.
				data = ::mxx2::scatterv(partition.begin(), partition.end(), send_counts, group, 0);

			}


			// scatter the ranges
			range = ::mxx2::scatter(ranges, group, 0);

			// now create the L2Blocks from the data  (reuse block)
			block.assign(&(data[0]), &(data[0]) + range.size(), range);

			TIMER_END(file, "open", data.size());



			//== reserve
			TIMER_START(file);
			// broadcast the estimated size
			mxx::datatype<size_t> size_dt;
			MPI_Bcast(&est_size, 1, size_dt.type(), 0, group);
			result.reserve(est_size);
			TIMER_END(file, "reserve", est_size);


			// == parse kmer/tuples iterator
			TIMER_START(file);


			read_block<SeqParser, KP>(block, result);

			TIMER_END(file, "read", result.size());
		}
		TIMER_REPORT_MPI(file, rank, _comm);



		INFOF("freeing group communicator");
		MPI_Comm_free(&group);
		INFOF("freeing group_leader communicator");
		if (group_leaders != MPI_COMM_NULL) MPI_Comm_free(&group_leaders);
		INFOF("DONE WITH communicator release");


		return result.size() - before;
	}


	std::vector<TupleType> find(std::vector<KmerType> &query) {
		return map.find(query);
	}

	std::vector< std::pair<KmerType, size_t> > count(std::vector<KmerType> &query) {
		return map.count(query);
	}

	void erase(std::vector<KmerType> &query) {
		map.erase(query);
	}


	template <typename Predicate>
	std::vector<TupleType> find_if(std::vector<KmerType> &query, Predicate const &pred) {
		return map.find_if(query, pred);
	}
	template <typename Predicate>
	std::vector<TupleType> find_if(Predicate const &pred) {
		return map.find_if(pred);
	}

	template <typename Predicate>
	std::vector<std::pair<KmerType, size_t> > count_if(std::vector<KmerType> &query, Predicate const &pred) {
		return map.count_if(query, pred);
	}

	template <typename Predicate>
	std::vector<std::pair<KmerType, size_t> > count_if(Predicate const &pred) {
		return map.count_if(pred);
	}


	template <typename Predicate>
	void erase_if(std::vector<KmerType> &query, Predicate const &pred) {
		map.erase_if(query, pred);
	}

	template <typename Predicate>
	void erase_if(Predicate const &pred) {
		map.erase_if(pred);
	}

	size_t local_size() {
		return map.local_size();
	}

	/**
	 * @tparam T 	input type may not be same as map's value types, so map need to provide overloads (and potentially with transform operators)
	 */
	 template <typename T>
	void insert(std::vector<T> &temp) {
		TIMER_INIT(build);

		TIMER_START(build);
		this->map.reserve(this->map.size() + temp.size());
		TIMER_END(build, "reserve", temp.size());


		// distribute
		TIMER_START(build);
		this->map.insert(temp);
		TIMER_END(build, "insert", this->map.local_size());


		TIMER_START(build);
		size_t m = this->map.update_multiplicity();
		TIMER_END(build, "multiplicity", m);

		TIMER_REPORT_MPI(build, this->commRank, this->comm);

	 }

	 // Note that KmerParserType may depend on knowing the Sequence Parser Type (e.g. provide quality score iterators)
	 //	Output type of KmerParserType may not match Map value type, in which case the map needs to do its own transform.
	 //     since Kmer template parameter is not explicitly known, we can't hard code the return types of KmerParserType.

	 /// convenience function for building index.
	 void build(const std::string & filename, MPI_Comm comm) {

	   ::std::vector<typename KmerParser::value_type> temp;

	   // file extension determines SeqParserType
	   std::string extension = ::bliss::utils::file::get_file_extension(filename);
	   std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
	   if (extension.compare("fastq") == 0) {
	     // default to including quality score iterators.
	     this->read_file<::bliss::io::FASTQParser<true> >(filename, temp, comm);
	   } else {
	     throw std::invalid_argument("input filename extension is not supported.");
	   }

		 this->insert(temp);
	 }


	 void build_with_mpi_subcomm(const std::string & filename, MPI_Comm comm) {
		 ::std::vector<typename KmerParser::value_type> temp;

	    // file extension determines SeqParserType
	     std::string extension = ::bliss::utils::file::get_file_extension(filename);
	     std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
	     if (extension.compare("fastq") == 0) {
	       // default to including quality score iterators.
	       this->read_file_mpi_subcomm<::bliss::io::FASTQParser<true> >(filename, temp, comm);
	     } else {
	       throw std::invalid_argument("input filename extension is not supported.");
	     }

		 this->insert(temp);
	 }

};



template <typename KmerType>
struct KmerParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = KmerType;


  /**
   * @brief generate kmers from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {

		static_assert(std::is_same<KmerType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
						"output type and output container value type are not the same");

		using Alphabet = typename KmerType::KmerAlphabet;

		/// converter from ascii to alphabet values
		using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

		/// kmer generation iterator
		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIterType>::value_type,
				KmerType>::value,
				"input iterator and output iterator's value types differ");

		// then compute and store into index (this will generate kmers and insert into index)
		if (read.seqBegin == read.seqEnd) return output_iter;

		//== set up the kmer generating iterators.
		KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
		KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);

		return ::std::copy(start, end, output_iter);

	}
};



template <typename TupleType>
struct KmerPositionTupleParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = TupleType;


  /**
   * @brief generate kmer-position pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {

		static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
						"output type and output container value type are not the same");
		static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

		using KmerType = typename std::tuple_element<0, TupleType>::type;
		using Alphabet = typename KmerType::KmerAlphabet;

		using IdType = typename std::tuple_element<1, TupleType>::type;
		static_assert(::std::is_same<typename SeqType::IdType, IdType>::value, "position type does not match for input and output iterators" );

		/// converter from ascii to alphabet values
		using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

		/// kmer generation iterator
		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
		/// kmer position iterator type
		using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

		/// combine kmer iterator and position iterator to create an index iterator type.
		using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, IdIterType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
				TupleType>::value,
				"input iterator and output iterator's value types differ");


		// then compute and store into index (this will generate kmers and insert into index)
		if (read.seqBegin == read.seqEnd) return output_iter;

		//== set up the kmer generating iterators.
		KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
		KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);

		//== set up the position iterators
		IdIterType id_start(read.id);
		IdIterType id_end(read.id);

		// ==== set up the zip iterators
		KmerIndexIterType index_start(start, id_start);
		KmerIndexIterType index_end(end, id_end);

		return ::std::copy(index_start, index_end, output_iter);
	}

};



template <typename TupleType>
struct KmerPositionQualityTupleParser {

    /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = TupleType;

  /**
   * @brief generate kmer-position-quality pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
   * @param read          sequence object, which has pointers to the raw byte array.
   * @param output_iter   output iterator pointing to insertion point for underlying container.
   * @return new position for output_iter
   * @tparam SeqType      type of sequence.  inferred.
   * @tparam OutputIt     output iterator type, inferred.
   */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {


		static_assert(SeqType::has_quality::value, "Sequence Parser needs to support quality scores");

		static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
				"output type and output container value type are not the same");

		static_assert(::std::tuple_size<TupleType>::value == 2, "kmer-pos-qual index data type should be a pair");

		using KmerType = typename std::tuple_element<0, TupleType>::type;
		using Alphabet = typename KmerType::KmerAlphabet;

		using KmerInfoType = typename std::tuple_element<1, TupleType>::type;
		static_assert(::std::tuple_size<KmerInfoType>::value == 2, "pos-qual index data type should be a pair");

		using QualType = typename std::tuple_element<1, KmerInfoType>::type;
		using IdType = typename std::tuple_element<0, KmerInfoType>::type;
		static_assert(::std::is_same<typename SeqType::IdType, IdType>::value, "position type does not match for input and output iterators" );


		/// converter from ascii to alphabet values
		using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

		/// kmer generation iterator
		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

		/// kmer position iterator type
		using IdIterType = bliss::iterator::SequenceIdIterator<IdType>;

		using QualIterType = bliss::index::QualityScoreGenerationIterator<typename SeqType::IteratorType, 21, bliss::index::Illumina18QualityScoreCodec<QualType> >;

		/// combine kmer iterator and position iterator to create an index iterator type.
		using KmerInfoIterType = bliss::iterator::ZipIterator<IdIterType, QualIterType>;
		static_assert(std::is_same<typename std::iterator_traits<KmerInfoIterType>::value_type,
				KmerInfoType>::value,
				"kmer info input iterator and output iterator's value types differ");


		using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, KmerInfoIterType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
				TupleType>::value,
				"input iterator and output iterator's value types differ");


		// then compute and store into index (this will generate kmers and insert into index)
		if (read.seqBegin == read.seqEnd || read.qualBegin == read.qualEnd) return output_iter;

		//== set up the kmer generating iterators.
		KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
		KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);

		//== set up the position iterators
		IdIterType id_start(read.id);
		IdIterType id_end(read.id);

		QualIterType qual_start(read.qualBegin);
		QualIterType qual_end(read.qualEnd);

		KmerInfoIterType info_start(id_start, qual_start);
		KmerInfoIterType info_end(id_end, qual_end);


		// ==== set up the zip iterators
		KmerIndexIterType index_start(start, info_start);
		KmerIndexIterType index_end(end, info_end);


		return ::std::copy(index_start, index_end, output_iter);
	}
};




template <typename TupleType>
struct KmerCountTupleParser {

  /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.
	using value_type = TupleType;

	/**
	 * @brief generate kmer-count pairs from 1 sequence.  result inserted into output_iter, which may be preallocated.
	 * @param read          sequence object, which has pointers to the raw byte array.
	 * @param output_iter   output iterator pointing to insertion point for underlying container.
	 * @return new position for output_iter
	 * @tparam SeqType      type of sequence.  inferred.
	 * @tparam OutputIt     output iterator type, inferred.
	 */
	template <typename SeqType, typename OutputIt>
	OutputIt operator()(SeqType & read, OutputIt output_iter) {

		static_assert(std::is_same<TupleType, typename ::std::iterator_traits<OutputIt>::value_type>::value,
				"output type and output container value type are not the same");
		static_assert(::std::tuple_size<TupleType>::value == 2, "count data type should be a pair");

		using KmerType = typename std::tuple_element<0, TupleType>::type;
		using CountType = typename std::tuple_element<1, TupleType>::type;

		using Alphabet = typename KmerType::KmerAlphabet;

		/// converter from ascii to alphabet values
		using BaseCharIterator = bliss::iterator::transform_iterator<typename SeqType::IteratorType, bliss::common::ASCII2<Alphabet> >;

		/// kmer generation iterator
		using KmerIterType = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;
		using CountIterType = bliss::iterator::ConstantIterator<CountType>;

		using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIterType, CountIterType>;

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type, std::pair<KmerType, CountType> >::value,
				"count zip iterator not producing the right type.");

		static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
				TupleType>::value,
				"count: input iterator and output iterator's value types differ");

		CountIterType count_start(1);

		// then compute and store into index (this will generate kmers and insert into index)
		if (read.seqBegin == read.seqEnd) return output_iter;

		//== set up the kmer generating iterators.
		KmerIterType start(BaseCharIterator(read.seqBegin, bliss::common::ASCII2<Alphabet>()), true);
		KmerIterType end(BaseCharIterator(read.seqEnd, bliss::common::ASCII2<Alphabet>()), false);

		KmerIndexIterType istart(start, count_start);
		KmerIndexIterType iend(end, count_start);

		return ::std::copy(istart, iend, output_iter);
	}

};

// version that computes the
//template <typename MapType>
//using KmerIndex = Index<MapType, KmerParser<typename MapType::key_type> >;

template <typename MapType>
using PositionIndex = Index<MapType, KmerPositionTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;

template <typename MapType>
using PositionQualityIndex = Index<MapType, KmerPositionQualityTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;

template <typename MapType>
using CountIndex = Index<MapType, KmerCountTupleParser<std::pair<typename MapType::key_type, typename MapType::mapped_type> > >;


} /* namespace kmer */

} /* namespace index */
} /* namespace bliss */

#endif /* KMERINDEX2_HPP_ */
