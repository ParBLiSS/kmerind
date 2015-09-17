/*
 * de_bruijn_construct_engine.hpp
 *
 *  Created on: Aug 17, 2015
 *      Author: yongchao
 */

#ifndef DE_BRUIJN_CONSTRUCT_ENGINE_HPP_
#define DE_BRUIJN_CONSTRUCT_ENGINE_HPP_


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
#include "io/fasta_loader.hpp"
//#include "io/fasta_iterator.hpp"

#include "utils/logging.h"
#include "common/alphabets.hpp"
#include "common/kmer.hpp"
#include "common/base_types.hpp"
#include "common/sequence.hpp"
#include "utils/kmer_utils.hpp"

#include "io/mxx_support.hpp"
#include "containers/distributed_unordered_map.hpp"
#include "containers/distributed_sorted_map.hpp"
// way too slow.  also not updated. #include <containers/distributed_map.hpp>

#include "io/sequence_iterator.hpp"
#include "io/sequence_id_iterator.hpp"
#include "iterators/transform_iterator.hpp"
#include "common/kmer_iterators.hpp"
#include "iterators/zip_iterator.hpp"
#include "iterators/constant_iterator.hpp"
#include "index/quality_score_iterator.hpp"
#include "wip/kmer_hash.hpp"
#include "wip/de_bruijn_node_trait.hpp"
#include "iterators/edge_iterator.hpp"
#include "wip/kmer_index.hpp"

#include "utils/timer.hpp"

namespace bliss
{
  namespace de_bruijn
  {
//    	/*generic class to create a directed de Bruijn graph*/
//      template <typename MapType, typename Parser, typename EdgeType>
//      class bidirected{
//        protected:
//
//          MapType map;
//
//          MPI_Comm comm;
//          int commSize;
//          int commRank;
//
//        public:
//          using KmerType = typename MapType::key_type;
//          using NodeTraitType = typename MapType::mapped_type;
//          using TupleType =      std::pair<KmerType, EdgeType>;
//
//          using Alphabet = typename KmerType::KmerAlphabet;
//
//          bidirected(MPI_Comm _comm, int _comm_size) : map(_comm, _comm_size), comm(_comm) {
//            MPI_Comm_size(_comm, &commSize);
//            MPI_Comm_rank(_comm, &commRank);
//          }
//
//          virtual ~bidirected() {/*do nothing*/};
//
//          MapType & get_map() {
//            return map;
//          }
//
//          /**
//           * @brief  generate kmers or kmer tuples for 1 block of raw data.
//           * @note   requires that SeqParser be passed in and operates on the Block's Iterators.
//           *          Mostly, this is because we need to broadcast the state of SeqParser to procs on the same node (L1 seq info) and recreate on child procs a new SeqParser (for L2)
//           * @tparam SeqParser    parser type for extracting sequences.  supports FASTQ and FASTA.  template template parameter, param is iterator
//           * @tparam KP           parser type for generating Kmer.  supports kmer, kmer+pos, kmer+count, kmer+pos/qual.
//           * @tparam BlockType    input partition type, supports in memory (vector) vs memmapped.
//           * @param partition
//           * @param result        output vector.  should be pre allocated.
//           */
//          template <typename KP, template <typename> class SeqParser, typename BlockType>
//          static size_t read_block(BlockType & partition, SeqParser<typename BlockType::iterator> const &seq_parser, std::vector<typename KP::value_type>& result) {
//
//            // from FileLoader type, get the block iter type and range type
//            using BlockIterType = typename BlockType::iterator;
//
//            using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, SeqParser >;
//            using SeqType = typename ::std::iterator_traits<SeqIterType>::value_type;
//
//            //== sequence parser type
//            KP kmer_parser;
//
//            //== process the chunk of data
//            SeqType read;
//
//            //==  and wrap the chunk inside an iterator that emits Reads.
//            SeqIterType seqs_start(seq_parser, partition.begin(), partition.end(), partition.getRange().start);
//            SeqIterType seqs_end(partition.end());
//
//            ::fsc::back_emplace_iterator<std::vector<typename KP::value_type> > emplace_iter(result);
//
//            size_t before = result.size();
//
//            //== loop over the reads
//            for (; seqs_start != seqs_end; ++seqs_start)
//            {
//        //      std::cout << "** seq: " << (*seqs_start).id.id << ", ";
//        //      ostream_iterator<typename std::iterator_traits<typename SeqType::IteratorType>::value_type> osi(std::cout);
//        //      std::copy((*seqs_start).seqBegin, (*seqs_start).seqEnd, osi);
//        //      std::cout << std::endl;
//
//              emplace_iter = kmer_parser(*seqs_start, emplace_iter);
//
//        //      std::cout << "Last: pos - kmer " << result.back() << std::endl;
//
//            }
//
//            return result.size() - before;
//          }
//
//
//
//          /*load reads from the database*/
//          template <template <typename> class SeqParser, typename KP = Parser>
//          static size_t read_file(const std::string & filename, std::vector<typename KP::value_type> & result, MPI_Comm _comm) {
//
//            int p, rank;
//            MPI_Comm_size(_comm, &p);
//            MPI_Comm_rank(_comm, &rank);
//
//
//            DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//            using FileLoaderType = bliss::io::FileLoader<CharType, 0, SeqParser, false, false>; // raw data type :  use CharType
//
//            //====  now process the file, one L1 block (block partition by MPI Rank) at a time
//            DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//            size_t before = result.size();
//            DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//            TIMER_INIT(file);
//            {  // ensure that fileloader is closed at the end.
//            	DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//              TIMER_START(file);
//              //==== create file Loader
//              FileLoaderType loader(filename, _comm, 1, sysconf(_SC_PAGE_SIZE));  // this handle is alive through the entire building process.
//              typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();
//              TIMER_END(file, "open", partition.getRange().size());
//              DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//
//
//              //== reserve
//              DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//              TIMER_START(file);
//              // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
//              // index reserve internally sends a message to itself.
//
//              // call after getting first L1Block to ensure that file is loaded.
//              DEBUG("P: " << p);
//              size_t est_size = (loader.getKmerCountEstimate(KmerType::size) + p - 1) / p;
//
//              DEBUG("est_size: " << est_size);
//              result.reserve(est_size);
//              TIMER_END(file, "reserve", est_size);
//              DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//
//              TIMER_START(file);
//              auto l1parser = loader.getSeqParser();
//              l1parser.init_for_iterator(partition.begin(), loader.getFileRange(), partition.getRange(), partition.getRange(), _comm);
//              TIMER_END(file, "mark_seqs", est_size);
//
//              TIMER_START(file);
//              //=== copy into array
//              while (partition.getRange().size() > 0) {
//                read_block<KP>(partition, l1parser, result);
//
//                partition = loader.getNextL1Block();
//              }
//              TIMER_END(file, "read", result.size());
//              DEBUG("function: " << __FUNCTION__ << " line: " << __LINE__);
//            }
//
//            TIMER_REPORT_MPI(file, rank, _comm);
//            return result.size() - before;
//          }
//
//
//          /*load reads from the database*/
//          template <template <typename> class SeqParser, typename KP = Parser>
//          static size_t read_file_mpi_subcomm(const std::string & filename, std::vector<typename KP::value_type>& result, MPI_Comm _comm) {
//
//            int p, rank;
//            MPI_Comm_size(_comm, &p);
//            MPI_Comm_rank(_comm, &rank);
//
//            // split the communcator so 1 proc from each host does the read, then redistribute.
//            MPI_Comm group_leaders;
//            MPI_Comm group;
//            ::std::tie(group_leaders, group) = ::mxx2::split_communicator_by_host(_comm);
//
//            int g_size, g_rank;
//            MPI_Comm_rank(group, &g_rank);
//            MPI_Comm_size(group, &g_size);
//
//            int gl_size = MPI_UNDEFINED;
//            int gl_rank = MPI_UNDEFINED;
//            if (group_leaders != MPI_COMM_NULL) {
//              MPI_Comm_rank(group_leaders, &gl_rank);
//              MPI_Comm_size(group_leaders, &gl_size);
//            }
//
//            // raw data type :  use CharType.   block partition at L1 and L2.  no buffering at all, since we will be copying data to the group members anyways.
//            using FileLoaderType = bliss::io::FileLoader<CharType, 0, SeqParser, false, false,
//                bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >,  bliss::partition::BlockPartitioner<bliss::partition::range<size_t> >>;
//
//            //====  now process the file, one L1 block (block partition by MPI Rank) at a time
//            size_t before = result.size();
//
//            typename FileLoaderType::RangeType file_range;
//
//            TIMER_INIT(file);
//            {
//              typename FileLoaderType::L2BlockType block;
//              typename FileLoaderType::RangeType range;
//              ::std::vector<CharType> data;
//
//              ::std::vector<typename FileLoaderType::RangeType> ranges;
//              ::std::vector<size_t> send_counts;
//              typename FileLoaderType::L1BlockType partition;
//
//              // first load the file using the group loader's communicator.
//              TIMER_START(file);
//              size_t est_size = 1;
//              if (group_leaders != MPI_COMM_NULL) {  // ensure file loader is closed properly.
//                //==== create file Loader. this handle is alive through the entire building process.
//                FileLoaderType loader(group_leaders, filename, g_size);  // for member of group_leaders, each create g_size L2blocks.
//
//                // modifying the local index directly here causes a thread safety issue, since callback thread is already running.
//
//                partition = loader.getNextL1Block();
//
//                // call after getting first L1Block to ensure that file is loaded.
//                est_size = (loader.getKmerCountEstimate(KmerType::size) + p - 1) / p;
//
//                //====  now compute the send counts and ranges to be scattered.
//                ranges.resize(g_size);
//                send_counts.resize(g_size);
//                for (size_t i = 0; i < g_size; ++i) {
//                  block = loader.getNextL2Block(i);
//
//                  ranges[i] = block.getRange();
//                  send_counts[i] = ranges[i].size();
//                }
//
//                // scatter the data .  this call here relies on FileLoader still having the memory mapped.
//                data = ::mxx2::scatterv(partition.begin(), partition.end(), send_counts, group, 0);
//
//                // send the file range.
//                file_range = loader.getFileRange();
//
//
//              } else {
//                // replicated here because the root's copy needs to be inside the if clause - requires FileLoader to be open for the sender.
//
//                // scatter the data .  this call here relies on FileLoader still having the memory mapped.
//                data = ::mxx2::scatterv(partition.begin(), partition.end(), send_counts, group, 0);
//
//              }
//
//
//              // scatter the ranges
//              mxx::datatype<typename FileLoaderType::RangeType> range_dt;
//              MPI_Bcast(&file_range, 1, range_dt.type(), 0, group);
//
//              range = ::mxx2::scatter(ranges, group, 0);
//
//              // now create the L2Blocks from the data  (reuse block)
//              block.assign(&(data[0]), &(data[0]) + range.size(), range);
//
//              TIMER_END(file, "open", data.size());
//
//              // not reusing the SeqParser in loader.  instead, reinitializing one.
//              TIMER_START(file);
//              SeqParser<typename FileLoaderType::L2BlockType::iterator> l2parser;
//              l2parser.init_for_iterator(block.begin(), file_range, range, range, _comm);
//              TIMER_END(file, "mark_seqs", est_size);
//
//              //== reserve
//              TIMER_START(file);
//              // broadcast the estimated size
//              mxx::datatype<size_t> size_dt;
//              MPI_Bcast(&est_size, 1, size_dt.type(), 0, group);
//              result.reserve(est_size);
//              TIMER_END(file, "reserve", est_size);
//
//
//              // == parse kmer/tuples iterator
//              TIMER_START(file);
//
//              read_block<KP>(block, l2parser, result);
//
//              TIMER_END(file, "read", result.size());
//            }
//            TIMER_REPORT_MPI(file, rank, _comm);
//
//
//
//            INFOF("freeing group communicator");
//            MPI_Comm_free(&group);
//            INFOF("freeing group_leader communicator");
//            if (group_leaders != MPI_COMM_NULL) MPI_Comm_free(&group_leaders);
//            INFOF("DONE WITH communicator release");
//
//
//            return result.size() - before;
//          }
//
//
//	std::vector<TupleType> find(std::vector<KmerType> &query) {
//            return map.find(query);
//          }
//
//	std::vector< std::pair<KmerType, size_t> > count(std::vector<KmerType> &query) {
//		return map.count(query);
//	}
//
//          void erase(std::vector<KmerType> &query) {
//            map.erase(query);
//          }
//
//          template <typename Predicate>
//	std::vector<TupleType> find_if(std::vector<KmerType> &query, Predicate const &pred) {
//            return map.find_if(query, pred);
//          }
//          template <typename Predicate>
//	std::vector<TupleType> find_if(Predicate const &pred) {
//            return map.find_if(pred);
//          }
//
//	template <typename Predicate>
//	std::vector<std::pair<KmerType, size_t> > count_if(std::vector<KmerType> &query, Predicate const &pred) {
//		return map.count_if(query, pred);
//	}
//
//	template <typename Predicate>
//	std::vector<std::pair<KmerType, size_t> > count_if(Predicate const &pred) {
//		return map.count_if(pred);
//	}
//
//
//          template <typename Predicate>
//          void erase_if(std::vector<KmerType> &query, Predicate const &pred) {
//            map.erase_if(query, pred);
//          }
//
//          template <typename Predicate>
//          void erase_if(Predicate const &pred) {
//            map.erase_if(pred);
//          }
//
//          size_t local_size() {
//            return map.local_size();
//          }
//
//
//          /*construct graph nodes from k-mers*/
//          template <typename T>
//          void insert(std::vector<T> &temp) {
//            TIMER_INIT(build);
//
//            TIMER_START(build);
//            this->map.reserve(this->map.size() + temp.size());
//            TIMER_END(build, "reserve", temp.size());
//
//
//            // distribute
//            TIMER_START(build);
//            this->map.insert(temp);
//            TIMER_END(build, "insert", this->map.local_size());
//
//
//            TIMER_START(build);
//            size_t m = this->map.update_multiplicity();
//            TIMER_END(build, "multiplicity", m);
//
//            TIMER_REPORT_MPI(build, this->commRank, this->comm);
//
//          }
//
//          /// convenience function for building index.
//          template <template <typename> class SeqParser>
//          void build(const std::string & filename, MPI_Comm comm) {
//
//            // file extension determines SeqParserType
//            std::string extension = ::bliss::utils::file::get_file_extension(filename);
//            std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
//            if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0)) {
//              throw std::invalid_argument("input filename extension is not supported.");
//            }
//
//            // check to make sure that the file parser will work
//            if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
//              throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
//            } else if ((extension.compare("fasta") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser2<char*> >::value)) {
//              throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
//            }
//
//            // proceed
//            ::std::vector<typename Parser::value_type> temp;
//            this->read_file<SeqParser >(filename, temp, comm);
//
//
//       //        // dump the generated kmers to see if they look okay.
//       //         std::stringstream ss;
//       //         ss << "test." << commRank << ".log";
//       //         std::ofstream ofs(ss.str());
//       //         for (int i = 0; i < temp.size(); ++i) {
//       //          ofs << "item: " << i << " value: " << temp[i] << std::endl;
//       //         }
//       //         ofs.close();
//
//            DEBUG("Last: pos - kmer " << temp.back());
//            this->insert(temp);
//
//          }
//
//          template <template <typename> class SeqParser>
//          void build_with_mpi_subcomm(const std::string & filename, MPI_Comm comm) {
//              // file extension determines SeqParserType
//              std::string extension = ::bliss::utils::file::get_file_extension(filename);
//              std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
//              if ((extension.compare("fastq") != 0) && (extension.compare("fasta") != 0)) {
//                throw std::invalid_argument("input filename extension is not supported.");
//              }
//
//              // check to make sure that the file parser will work
//              if ((extension.compare("fastq") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTQParser<char*> >::value)) {
//                throw std::invalid_argument("Specified File Parser template parameter does not support files with fastq extension.");
//              } else if ((extension.compare("fasta") == 0) && (!std::is_same<SeqParser<char*>, ::bliss::io::FASTAParser2<char*> >::value)) {
//                throw std::invalid_argument("Specified File Parser template parameter does not support files with fasta extension.");
//              }
//
//              // proceed
//              ::std::vector<typename Parser::value_type> temp;
//              this->read_file_mpi_subcomm<SeqParser >(filename, temp, comm);
//
//              DEBUG("Last: pos - kmer " << temp.back());
//              this->insert(temp);
//
//          }
//
//      };

  /*generate de Bruijn graph nodes and edges from the reads*/
    /**
     * @tparam TupleType  value type of outputIt, or the generated type.  supplied so that we can use template template param with Index.
     */
   template <typename KmerType, typename EdgeEncoder = bliss::common::DNA16>
	 struct de_bruijn_parser {

      /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUTIT.
      using edge_type = typename std::conditional<std::is_same<EdgeEncoder, bliss::common::ASCII>::value,
                                                  uint16_t, uint8_t>::type;
      using value_type = std::pair<KmerType, edge_type>;

      using Alphabet = typename KmerType::KmerAlphabet;

		 template <typename SeqType, typename OutputIt>
		 OutputIt operator()(SeqType & read, OutputIt output_iter) {
		   static_assert(std::is_same<value_type, typename ::std::iterator_traits<OutputIt>::value_type>::value,
		                  "output type and output container value type are not the same");

	      // filter out EOL characters
	       using CharIter = bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>;
	       // converter from ascii to alphabet values
	       using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
	       // kmer generation iterator
	       using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;


		    // edge iterator type
//		   using EdgeType = typename ::std::tuple_element<1, TupleType>::type;
//		   static_assert((std::is_same<EdgeEncoder, bliss::common::ASCII>::value && std::is_same<EdgeType, uint16_t>::value) ||
//		                 std::is_same<EdgeType, uint8_t>::value, "Edge Type is not compatible with EdgeEncoder type");

       using EdgeIterType = bliss::iterator::edge_iterator<CharIter, EdgeEncoder>;


		   // combine kmer iterator and position iterator to create an index iterator type.
		   using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, EdgeIterType>;

		    static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
		                  value_type>::value,
		        "Generating iterator and output iterator's value types differ");

		    // then compute and store into index (this will generate kmers and insert into index)
		     if (read.seqBegin == read.seqEnd) return output_iter;

		     //== set up the kmer generating iterators.
		     bliss::index::kmer::NotEOL neol;
		     KmerIter start(BaseCharIterator(CharIter(neol, read.seqBegin, read.seqEnd), bliss::common::ASCII2<Alphabet>()), true);
		     KmerIter end(BaseCharIterator(CharIter(neol, read.seqEnd), bliss::common::ASCII2<Alphabet>()), false);


		     // set up edge iterator
	       EdgeIterType edge_start(CharIter(neol, read.seqBegin, read.seqEnd), CharIter(neol, read.seqEnd), KmerType::size);
	       EdgeIterType edge_end (CharIter(neol, read.seqEnd));


	       // ==== set up the zip iterators
	       KmerIndexIterType node_start (start, edge_start);
	       KmerIndexIterType node_end(end, edge_end);

	       return ::std::copy(node_start, node_end, output_iter);

		 }
	 };

	 /*generate de Brujin graph nodes and edges, which each node associated with base quality scores*/
    template <typename KmerType, typename QualType=double, typename EdgeEncoder = bliss::common::DNA16, template<typename> class QualityEncoder = bliss::index::Illumina18QualityScoreCodec>
	 struct de_bruijn_quality_parser {

        /// type of element generated by this parser.  since kmer itself is parameterized, this is not hard coded.  NOTE THAT THIS IS TYPE FOR THE OUTPUTIT.
        using edge_type = typename std::conditional<std::is_same<EdgeEncoder, bliss::common::ASCII>::value,
                                                    uint16_t, uint8_t>::type;
        using value_type = std::pair<KmerType, std::pair<edge_type, QualType> >;

        using Alphabet = typename KmerType::KmerAlphabet;

       template <typename SeqType, typename OutputIt>
       OutputIt operator()(SeqType & read, OutputIt output_iter) {

         static_assert(std::is_same<value_type, typename ::std::iterator_traits<OutputIt>::value_type>::value,
                        "output type and output container value type are not the same");


         // filter out EOL characters
          using CharIter = bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>;
          // converter from ascii to alphabet values
          using BaseCharIterator = bliss::iterator::transform_iterator<CharIter, bliss::common::ASCII2<Alphabet> >;
          // kmer generation iterator
          using KmerIter = bliss::common::KmerGenerationIterator<BaseCharIterator, KmerType>;

          // edge iterator type
//          using EdgeType = typename std::tuple_element<0, typename std::tuple_element<1, TupleType>::type >::type;
//          static_assert((std::is_same<EdgeEncoder, bliss::common::ASCII>::value && std::is_same<EdgeType, uint16_t>::value) ||
//                          std::is_same<EdgeType, uint8_t>::value, "Edge Type is not compatible with EdgeEncoder type");

          using EdgeIterType = bliss::iterator::edge_iterator<CharIter, EdgeEncoder>;

          // also remove eol from quality score
          using QualIterType =
              bliss::index::QualityScoreGenerationIterator<bliss::index::kmer::NonEOLIter<typename SeqType::IteratorType>, KmerType::size, QualityEncoder<QualType> >;

          // combine kmer iterator and position iterator to create an index iterator type.
          using KmerInfoIterType = bliss::iterator::ZipIterator<EdgeIterType, QualIterType>;

          using KmerIndexIterType = bliss::iterator::ZipIterator<KmerIter, KmerInfoIterType>;

          static_assert(std::is_same<typename std::iterator_traits<KmerIndexIterType>::value_type,
                                     value_type>::value,
              "generating iterator and output iterator's value types differ");


           // then compute and store into index (this will generate kmers and insert into index)
           if (read.seqBegin == read.seqEnd || read.qualBegin == read.qualEnd) return output_iter;

           //== set up the kmer generating iterators.
           bliss::index::kmer::NotEOL neol;
           KmerIter start(BaseCharIterator(CharIter(neol, read.seqBegin, read.seqEnd), bliss::common::ASCII2<Alphabet>()), true);
           KmerIter end(BaseCharIterator(CharIter(neol, read.seqEnd), bliss::common::ASCII2<Alphabet>()), false);


           // set up edge iterator
           EdgeIterType edge_start(CharIter(neol, read.seqBegin, read.seqEnd), CharIter(neol, read.seqEnd), KmerType::size);
           EdgeIterType edge_end (CharIter(neol, read.seqEnd));


           QualIterType qual_start(CharIter(neol, read.qualBegin, read.qualEnd));
            QualIterType qual_end(CharIter(neol, read.qualEnd));

            KmerInfoIterType info_start(edge_start, qual_start);
            KmerInfoIterType info_end(edge_end, qual_end);



           // ==== set up the zip iterators
            KmerIndexIterType node_start(start, info_start);
            KmerIndexIterType node_end(end, info_end);

           return ::std::copy(node_start, node_end, output_iter);

		 }
	 };

	template <typename MapType, typename EdgeEncoder = bliss::common::DNA16>
	using de_bruijn_engine = ::bliss::index::kmer::Index<MapType, de_bruijn_parser<typename MapType::key_type, EdgeEncoder > >;

	template <typename MapType, typename EdgeEncoder = bliss::common::DNA16, template<typename> class QualityEncoder = bliss::index::Illumina18QualityScoreCodec >
  using de_bruijn_quality_engine = ::bliss::index::kmer::Index<MapType, de_bruijn_quality_parser<typename MapType::key_type, typename std::tuple_element<1, typename MapType::mapped_type>::type, EdgeEncoder, QualityEncoder > >;

  } /* namespace de_bruijn */
} /* namespace bliss */

#endif /* DE_BRUIJN_CONSTRUCT_ENGINE_HPP_ */
