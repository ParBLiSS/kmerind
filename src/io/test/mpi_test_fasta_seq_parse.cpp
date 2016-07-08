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
 * fastaloader_test.cpp
 *
 *  Created on: July 4, 2014
 *      Author: Tony Pan <tpan7@gatech.edu>
 */


#include "bliss-config.hpp"    // for location of data.

#if defined(USE_MPI)
#include "mxx/env.hpp"
#include <mxx/reduction.hpp>
#endif

// include google test
#include <gtest/gtest.h>
#include <cstdint> // for uint64_t, etc.
#include <string>
#include <limits>

#include "common/sequence.hpp"
#include "io/sequence_iterator.hpp"

#include "index/kmer_index.hpp"

#include "io/test/file_load_test_fixtures.hpp"

#include "io/fasta_loader.hpp"
#include "io/file.hpp"




class FASTAParseTest : public KmerReaderTest
{
	// can't compare the sequences directly since offsets may be different.
protected:

	template <typename ITER>
	using ParserType = bliss::io::FASTAParser<ITER>;

	static constexpr size_t kmer_size = 35;
	using KmerType = ::bliss::common::Kmer<kmer_size, bliss::common::DNA5, uint64_t>;

	template <typename file_loader>
	void parse_seq(file_loader & fobj ) {
		  // get fileName

			// get fileName
			bliss::io::file_data fdata = fobj.read_file();

			ASSERT_TRUE(fobj.size() > 0);

			ASSERT_EQ(fdata.getRange().size(), fobj.size());
			ASSERT_EQ(fdata.getRange().end, fobj.size());

			ASSERT_EQ(fdata.in_mem_range_bytes.size(), fobj.size());
			ASSERT_EQ(fdata.in_mem_range_bytes.end, fobj.size());

			ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
			ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

		  ASSERT_TRUE((fdata.data[0] == '>') || (fdata.data[0] == ';'));

			this->seqCount = 0;
			this->kmerCount = 0;

			// from FileLoader type, get the block iter type and range type
			using BlockIterType = typename bliss::io::file_data::const_iterator;
			using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

			ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

			size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes);

			  //==  and wrap the chunk inside an iterator that emits Reads.
			  SeqIterType seqs_start(l1parser, fdata.begin(),
					  fdata.in_mem_end(), offset);
			  SeqIterType seqs_end(fdata.in_mem_end());

		      bliss::index::kmer::KmerParser<KmerType > kmer_parser;
		      std::vector<KmerType> result;

		      bool localcomp = true, comp = true;

		      std::vector<unsigned char> gold;

		      ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);


			  //== loop over the reads
			  for (; seqs_start != seqs_end; ++seqs_start)
			  {
//				  std::cout << " seq: " << *(seqs_start) << std::endl;


//				if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) {
					++(this->seqCount);

			          result.clear();

			          // generate kmers and save them
			          kmer_parser(*seqs_start, emplace_iter);

			          gold.clear();
			          gold.resize(seqs_start->seq_size(), 0);
			          this->readFilePOSIX(this->fileName, seqs_start->seq_global_offset(), seqs_start->seq_size(), gold.data());

			          // compare the results.
			          for (size_t i = 0; i < result.size(); ++i) {

			            localcomp = equal(bliss::utils::KmerUtils::toASCIIString(result[i]).begin(), &(gold[i]), kmer_size);

				          if (!localcomp) {
				        	  std::cout << "kmer: " << bliss::utils::KmerUtils::toASCIIString(result[i]) << std::endl;
				        	  std::cout << "gold: ";
				        	  for (size_t j = 0; j < kmer_size; ++j) {
				        		  std::cout << gold[i+j];
				        	  }
				        	  std::cout << std::endl;
				          }

				          comp &= localcomp;

			          }  // end kmers for



			          this->kmerCount += result.size();

			          ASSERT_TRUE(localcomp);

				}
		  //            std::cout << *seqs_start << ", ";
		  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
				//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
//			  }

	}


	template <typename file_loader>
	void parse_mpi(file_loader & fobj, size_t const & overlap, mxx::comm const & comm) {

		  // get this->fileName

		::bliss::io::file_data fdata = fobj.read_file();

		ASSERT_TRUE(fobj.size() > 0);

		// parent range should match.
		ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
		ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

		// get the ranges to makes sure they are in memory.
		ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

		ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);

//		std::cout << "rank " << comm.rank() << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
		// get all the sizes to make sure total is same as file size.
		size_t region_size = fdata.valid_range_bytes.size();
		region_size = ::mxx::allreduce(region_size, comm);
		ASSERT_EQ(region_size, fobj.size());

		// make sure the aggregate range is same as parent range.
		std::vector<size_t> begins = mxx::allgather(fdata.valid_range_bytes.start, comm);
		std::vector<size_t> ends = mxx::allgather(fdata.valid_range_bytes.end, comm);

		ASSERT_EQ(begins.front(), 0UL);
		ASSERT_EQ(ends.back(), fobj.size());

		int j = 0;
		for (int i = 1; i < comm.size(); ++i) {
//			if (comm.rank() == 0) std::cout << "rank " << comm.rank() << " end "  << (i-1) << ": " << ends[i-1] << " begin next " << i << ": " << begins[i] << " overlap = " << overlap << std::endl;
			if (begins[i] == ends[i]) continue;

			ASSERT_TRUE(ends[j] == (begins[i]) );  // valid range does not include overlap
			j = i;
		}
	//	  std::cout << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

		  this->seqCount = 0;
			this->kmerCount = 0;

		// from FileLoader type, get the block iter type and range type
		using BlockIterType = typename ::bliss::io::file_data::const_iterator;
		using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTAParser >;

		ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

		size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes, comm);

		  //==  and wrap the chunk inside an iterator that emits Reads.
		  SeqIterType seqs_start(l1parser, fdata.begin(),
				  fdata.in_mem_end(), offset);
		  SeqIterType seqs_end(fdata.in_mem_end());

	      bliss::index::kmer::KmerParser<KmerType > kmer_parser;
	      std::vector<KmerType> result;

	      bool localcomp = true, comp = true;

	      std::vector<unsigned char> gold;

	      ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);

		  //== loop over the reads
		  for (; seqs_start != seqs_end; ++seqs_start)
		  {
	          //std::cout << "rank " << comm.rank() << " seq: " << *(seqs_start) << std::endl;

	          // only count sequences that started in this partition and has non-zero size.
				if (((*seqs_start).record_size > 0) && ((*seqs_start).id.get_pos() >= fdata.valid_range_bytes.start)) ++(this->seqCount);

				// but generate kmers from all...
		          result.clear();

		          // generate kmers and save them
		          kmer_parser(*seqs_start, emplace_iter);
		          // std::cout << "rank " << comm.rank() << " seqs " << this->seqCount << " kmers: " << result.size() << std::endl;

		          gold.clear();
		          gold.resize(seqs_start->seq_size(), 0);
		          this->readFilePOSIX(this->fileName, seqs_start->seq_global_offset(), seqs_start->seq_size(), gold.data());

		          // compare the results.
		          for (size_t i = 0; i < result.size(); ++i) {

		            localcomp = equal(bliss::utils::KmerUtils::toASCIIString(result[i]).begin(), &(gold[i]), kmer_size);

			          if (!localcomp) {
			        	  std::cout << "kmer: " << bliss::utils::KmerUtils::toASCIIString(result[i]) << std::endl;
			        	  std::cout << "gold: ";
			        	  for (size_t j = 0; j < kmer_size; ++j) {
			        		  std::cout << gold[i+j];
			        	  }
			        	  std::cout << std::endl;
			          }

			          comp &= localcomp;
		          }  // end kmers for

//		          if (comm.rank() == 0) std::cout << "rank " << comm.rank() << " total count " << this->kmerCount << " + " << result.size() << std::endl;

		          this->kmerCount += result.size();

		          ASSERT_TRUE(comp);

			}
	  //            std::cout << *seqs_start << ", ";
	  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
			//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
//		  }
	//	  int rank;
	//	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//	  std::cout << "rank " << rank << " elem count " << this->seqCount << std::endl;

	//	  if (this->seqCount == 0) {
	//		  std::cout << " mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
	//	  }

	}

};


TEST_P(FASTAParseTest, parse_mmap)
{
#ifdef USE_MPI
	  ::mxx::comm comm;
	if (comm.rank() == 0) {
#endif

		bliss::io::mmap_file fobj(this->fileName);

		this->parse_seq(fobj);

#ifdef USE_MPI
	}
#endif

}

TEST_P(FASTAParseTest, parse_stdio)
{
#ifdef USE_MPI
	  ::mxx::comm comm;
	if (comm.rank() == 0) {
#endif

		bliss::io::stdio_file fobj(this->fileName);

		this->parse_seq(fobj);

#ifdef USE_MPI
	}
#endif

}

TEST_P(FASTAParseTest, parse_posix)
{
#ifdef USE_MPI
	  ::mxx::comm comm;
	if (comm.rank() == 0) {
#endif

    bliss::io::posix_file fobj(this->fileName);

    this->parse_seq(fobj);

#ifdef USE_MPI
  }
#endif

}

#ifdef USE_MPI
TEST_P(FASTAParseTest, parse_mmap_mpi)
{
	  ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, ParserType> fobj(this->fileName, kmer_size - 1, comm);

  this->parse_mpi(fobj, kmer_size - 1, comm);

  comm.barrier();
}

TEST_P(FASTAParseTest, parse_stdio_mpi)
{
	  ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, ParserType> fobj(this->fileName, kmer_size - 1, comm);

  this->parse_mpi(fobj, kmer_size - 1, comm);

  comm.barrier();
}

TEST_P(FASTAParseTest, parse_posix_mpi)
{
	  ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, ParserType> fobj(this->fileName, kmer_size - 1, comm);

  this->parse_mpi(fobj, kmer_size - 1, comm);

  comm.barrier();
}

TEST_P(FASTAParseTest, parse_mpiio_mpi)
{
	  ::mxx::comm comm;

	::bliss::io::parallel::mpiio_file<bliss::io::FASTAParser> fobj(this->fileName, kmer_size - 1, comm);

	  this->parse_mpi(fobj, kmer_size - 1, comm);

	  comm.barrier();
}
#endif




INSTANTIATE_TEST_CASE_P(Bliss, FASTAParseTest, ::testing::Values(
    TestFileInfo(243,  434,    13790, std::string("/test/data/natural.fasta")),
    TestFileInfo(250,  500,    14625, std::string("/test/data/natural.withN.fasta")),
    TestFileInfo(6,    246,    940, std::string("/test/data/test2.fasta")),
    TestFileInfo(4,    64,     512, std::string("/test/data/test.fasta")),
    TestFileInfo(5000, 335000, 1092580, std::string("/test/data/test.medium.fasta")),
    TestFileInfo(2,    31155,  31245, std::string("/test/data/test.unitiqs.fasta"))
));

//INSTANTIATE_TEST_CASE_P(Bliss, FASTAParseProcedureTest, ::testing::Values(
//		TestFileInfo(250, 14625, std::string("/test/data/natural.fasta"))
//		));

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

#if defined(USE_MPI)
  MPI_Init(&argc, &argv);
#endif

  result = RUN_ALL_TESTS();

#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  return result;
}

