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
 * fastqloader_test.cpp
 *
 *  Created on: Jun 3, 2014
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
#include <chrono>

#include "common/sequence.hpp"
#include "io/sequence_iterator.hpp"
#include "io/filtered_sequence_iterator.hpp"


#include "index/kmer_index.hpp"

#include "io/test/file_load_test_fixtures.hpp"

#include "io/fastq_loader.hpp"
#include "io/file.hpp"

#include "utils/benchmark_utils.hpp"




class FASTQParseTest : public KmerReaderTest
{
	// can't compare the sequences directly since offsets may be different.
protected:
	template <typename ITER>
	using ParserType = bliss::io::FASTQParser<ITER>;

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

		  ASSERT_TRUE(fdata.data[0] == '@');

			this->seqCount = 0;
			this->kmerCount = 0;

			// from FileLoader type, get the block iter type and range type
			using BlockIterType = typename ::bliss::io::file_data::const_iterator;
			using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

			ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

			size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes);

			  //==  and wrap the chunk inside an iterator that emits Reads.
			  SeqIterType seqs_start(l1parser, fdata.begin(),
					  fdata.in_mem_end(), offset);
			  SeqIterType seqs_end(fdata.in_mem_end());

		      bliss::index::kmer::KmerParser<KmerType > kmer_parser(fdata.valid_range_bytes);
		      std::vector<KmerType> result;

		      bool localcomp = true, comp = true;

		      std::vector<unsigned char> gold;

		      ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);


		      std::chrono::steady_clock::time_point t1, t2;
		      double duration = 0.0;
		      std::chrono::duration<double> time_span;

			  //== loop over the reads
			  for (; seqs_start != seqs_end; ++seqs_start)
			  {
				if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) {
					++(this->seqCount);

			          t1 = std::chrono::steady_clock::now();
			          result.clear();

			          // generate kmers and save them
			          kmer_parser(*seqs_start, emplace_iter);
			          t2 = std::chrono::steady_clock::now();
			          time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1));
			          duration += time_span.count();

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
			  }
			  std::cout << "parse time = " << duration << " s" << std::endl;

	}


	template <typename file_loader>
	void parse_mpi(file_loader & fobj, size_t const & overlap, mxx::comm const & comm) {

		  // get this->fileName

		::bliss::io::file_data fdata = fobj.read_file();

//		std::cout << "rank " << comm.rank() << fobj.get_class_name() << std::endl;

//		comm.barrier();
//        if (comm.rank() == 1) {
//        std::cout << "rank " << comm.rank() << "   parent " << fdata.parent_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   inmem " << fdata.in_mem_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   search " << fdata.valid_range_bytes << std::endl;
//
//			 std::cout << "rank " << comm.rank() << " ";
//			 for (auto ii = fdata.data.begin(); ii != (fdata.data.begin() + (fdata.in_mem_range_bytes.end - fdata.in_mem_range_bytes.start)); ++ii ) {
//				 std::cout << static_cast<unsigned char>(*ii);
//			 }
//			 std::cout << std::endl;
//        }
//        comm.barrier();

//		std::cout << " rank " << comm.rank() << " in mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

		ASSERT_TRUE(fobj.size() > 0);

		// parent range should match.
		ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
		ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

		// get the ranges to makes sure they are in memory.
		ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

		ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
		ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);

		if (fdata.valid_range_bytes.size() > 0) {
			if (fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] != '@') {
				std::cout << "ERROR: first char should be @, but is " << static_cast<unsigned char>(fdata.data[0]) << std::endl;
			}
			ASSERT_TRUE(fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] == '@');
		}

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
//			if (comm.rank() == 0) std::cout << "end "  << (i-1) << ": " << ends[i-1] << " begin next " << i << ": " << begins[i] << " overlap = " << overlap << std::endl;
			if (begins[i] == ends[i]) continue;

			ASSERT_TRUE(ends[j] == (begins[i]) );  // valid ranges do not include overlap.
			j = i;
		}
	//	  std::cout << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

		  this->seqCount = 0;
			this->kmerCount = 0;

		// from FileLoader type, get the block iter type and range type
		using BlockIterType = typename ::bliss::io::file_data::const_iterator;
		using SeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;

		ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

		size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes, comm);

		  //==  and wrap the chunk inside an iterator that emits Reads.
		  SeqIterType seqs_start(l1parser, fdata.begin(),
				  fdata.in_mem_end(), offset);
		  SeqIterType seqs_end(fdata.in_mem_end());

//          std::cout << "rank " << comm.rank() << " offset = " << offset << " sequence " << *seqs_start << std::endl;


	      bliss::index::kmer::KmerParser<KmerType > kmer_parser(fdata.valid_range_bytes);
	      std::vector<KmerType> result;

	      bool localcomp = true, comp = true;

	      std::vector<unsigned char> gold;

	      ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);

	      std::chrono::steady_clock::time_point t1, t2;
	      double duration = 0.0, gold_duration = 0.0, compare_duration = 0.0;
	      std::chrono::duration<double> time_span;

		  //== loop over the reads
		  for (; seqs_start != seqs_end; ++seqs_start)
		  {
			if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) {
				++(this->seqCount);

		          t1 = std::chrono::steady_clock::now();

		          result.clear();

//		          std::cout << "rank " << comm.rank() << " sequence " << *seqs_start << std::endl;

		          // generate kmers and save them
		          kmer_parser(*seqs_start, emplace_iter);
		          t2 = std::chrono::steady_clock::now();
		          time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1));
		          duration += time_span.count();

		          gold.clear();
		          gold.resize(seqs_start->seq_size(), 0);
		          this->readFilePOSIX(this->fileName, seqs_start->seq_global_offset(), seqs_start->seq_size(), gold.data());
	//	          std::cout << "rank " << comm.rank() << " entire gold: ";
	//	          for (size_t j = 0; j < seqs_start->seq_size(); ++j) {
	//	          			        		  std::cout << gold[j];
	//	          			        	  }
	//	          std::cout << std::endl;

		          t1 = std::chrono::steady_clock::now();
		          time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(t1 - t2));
		          gold_duration += time_span.count();


		          // compare the results.
		          for (size_t i = 0; i < result.size(); ++i) {

		            localcomp = equal(bliss::utils::KmerUtils::toASCIIString(result[i]).begin(), &(gold[i]), kmer_size);

			          if (!localcomp) {

			        	  std::cout << "rank " << comm.rank() << " kmer: " << bliss::utils::KmerUtils::toASCIIString(result[i]) << std::endl;
			        	  std::cout << "rank " << comm.rank() << " gold: ";
			        	  for (size_t j = 0; j < kmer_size; ++j) {
			        		  std::cout << gold[i+j];
			        	  }
			        	  std::cout << std::endl;
			          }

			          comp &= localcomp;
		          }  // end kmers for

		          this->kmerCount += result.size();
		          t2 = std::chrono::steady_clock::now();
		          time_span = (std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1));
		          compare_duration += time_span.count();


		          ASSERT_TRUE(comp);

			}
	  //            std::cout << *seqs_start << ", ";
	  //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
			//      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
		  }

		  std::cout << "rank " << comm.rank() << "parse time = " << duration << "s, gold parse time = " << gold_duration << "s, compare time = " << compare_duration << "s" << std::endl;
	//	  int rank;
	//	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//	  std::cout << "rank " << rank << " elem count " << this->seqCount << std::endl;

	//	  if (this->seqCount == 0) {
	//		  std::cout << " mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
	//	  }

	}

};


TEST_P(FASTQParseTest, parse_mmap)
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

TEST_P(FASTQParseTest, parse_stdio)
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

TEST_P(FASTQParseTest, parse_posix)
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
TEST_P(FASTQParseTest, parse_mmap_mpi)
{
	  ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(FASTQParseTest, parse_stdio_mpi)
{
	  ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(FASTQParseTest, parse_posix_mpi)
{
	  ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(FASTQParseTest, parse_mpiio_mpi)
{
	  ::mxx::comm comm;

	::bliss::io::parallel::mpiio_file<bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

	  this->parse_mpi(fobj, 0UL, comm);

	  comm.barrier();
}
#endif



// using 21-mer, first entry is number of records, second entry is file size.
INSTANTIATE_TEST_CASE_P(Bliss, FASTQParseTest, ::testing::Values(
    TestFileInfo(243,    434,     27580, std::string("/test/data/natural.fastq")),
    TestFileInfo(250,    500,     29250, std::string("/test/data/natural.withN.fastq")),
    TestFileInfo(7,      182,     939, std::string("/test/data/test.debruijn.small.fastq")),
    TestFileInfo(1,      26,      134, std::string("/test/data/test.debruijn.tiny.fastq")),
//    TestFileInfo(254562, 6618612, 34111308, std::string("/test/data/test.fastq")),
    TestFileInfo(140,    3640,    18761, std::string("/test/data/test.medium.fastq")),
    TestFileInfo(7,      182,     939, std::string("/test/data/test.small.fastq")),
    TestFileInfo(1,      16552,   33194, std::string("/test/data/test.unitiq1.fastq")),
    TestFileInfo(1,      3527,    7144, std::string("/test/data/test.unitiq1.short2.fastq")),
    TestFileInfo(1,      7367,    14824, std::string("/test/data/test.unitiq1.short.fastq")),
    TestFileInfo(1,      14603,   29296, std::string("/test/data/test.unitiq2.fastq")),
    TestFileInfo(2,      31155,   62490, std::string("/test/data/test.unitiqs.fastq"))
));




class NoNFASTQParseTest : public KmerReaderTest
{
  // can't compare the sequences directly since offsets may be different.
protected:
  template <typename ITER>
  using ParserType = bliss::io::FASTQParser<ITER>;

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

      ASSERT_TRUE(fdata.data[0] == '@');

      this->seqCount = 0;
      this->kmerCount = 0;

      // from FileLoader type, get the block iter type and range type
      using BlockIterType = typename ::bliss::io::file_data::const_iterator;
//    using OrigSeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;
//    using SeqIterType = ::bliss::iterator::filter_iterator<::bliss::io::SequenceNPredicate, OrigSeqIterType>;
      using SeqIterType = ::bliss::io::NFilterSequencesIterator<BlockIterType, bliss::io::FASTQParser >;

      ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

      size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes);

        //==  and wrap the chunk inside an iterator that emits Reads.
//        SeqIterType seqs_start(OrigSeqIterType(l1parser, fdata.begin(),
//            fdata.in_mem_end(), offset), OrigSeqIterType(fdata.in_mem_end()));
//        SeqIterType seqs_end(OrigSeqIterType(fdata.in_mem_end()));
        SeqIterType seqs_start(l1parser, fdata.begin(), fdata.in_mem_end(), offset);
        SeqIterType seqs_end(fdata.in_mem_end());

          bliss::index::kmer::KmerParser<KmerType > kmer_parser(fdata.valid_range_bytes);
          std::vector<KmerType> result;

          bool localcomp = true, comp = true;

          std::vector<unsigned char> gold;

          ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);


        //== loop over the reads
        for (; seqs_start != seqs_end; ++seqs_start)
        {
        if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) {
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
        }

  }


  template <typename file_loader>
  void parse_mpi(file_loader & fobj, size_t const & overlap, mxx::comm const & comm) {

      // get this->fileName

    ::bliss::io::file_data fdata = fobj.read_file();

//    std::cout << "rank " << comm.rank() << fobj.get_class_name() << std::endl;

//    comm.barrier();
//        if (comm.rank() == 1) {
//        std::cout << "rank " << comm.rank() << "   parent " << fdata.parent_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   inmem " << fdata.in_mem_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   search " << fdata.valid_range_bytes << std::endl;
//
//       std::cout << "rank " << comm.rank() << " ";
//       for (auto ii = fdata.data.begin(); ii != (fdata.data.begin() + (fdata.in_mem_range_bytes.end - fdata.in_mem_range_bytes.start)); ++ii ) {
//         std::cout << static_cast<unsigned char>(*ii);
//       }
//       std::cout << std::endl;
//        }
//        comm.barrier();

//    std::cout << " rank " << comm.rank() << " in mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

    ASSERT_TRUE(fobj.size() > 0);

    // parent range should match.
    ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
    ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

    // get the ranges to makes sure they are in memory.
    ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
    ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

    ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
    ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);

    if (fdata.valid_range_bytes.size() > 0) {
      if (fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] != '@') {
        std::cout << "ERROR: first char should be @, but is " << static_cast<unsigned char>(fdata.data[0]) << std::endl;
      }
      ASSERT_TRUE(fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] == '@');
    }

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
//      if (comm.rank() == 0) std::cout << "end "  << (i-1) << ": " << ends[i-1] << " begin next " << i << ": " << begins[i] << " overlap = " << overlap << std::endl;
      if (begins[i] == ends[i]) continue;

      ASSERT_TRUE(ends[j] == (begins[i]) );  // valid ranges do not include overlap.
      j = i;
    }
  //    std::cout << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

      this->seqCount = 0;
      this->kmerCount = 0;

    // from FileLoader type, get the block iter type and range type
    using BlockIterType = typename ::bliss::io::file_data::const_iterator;
//    using OrigSeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;
//    using SeqIterType = ::bliss::iterator::filter_iterator<::bliss::io::SequenceNPredicate, OrigSeqIterType>;
    using SeqIterType = ::bliss::io::NFilterSequencesIterator<BlockIterType, bliss::io::FASTQParser >;

    ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

    size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes, comm);

      //==  and wrap the chunk inside an iterator that emits Reads.
//      SeqIterType seqs_start(OrigSeqIterType(l1parser, fdata.begin(),
//          fdata.in_mem_end(), offset), OrigSeqIterType(fdata.in_mem_end()));
//      SeqIterType seqs_end(OrigSeqIterType(fdata.in_mem_end()));
    SeqIterType seqs_start(l1parser, fdata.begin(), fdata.in_mem_end(), offset);
    SeqIterType seqs_end(fdata.in_mem_end());


//          std::cout << "rank " << comm.rank() << " offset = " << offset << " sequence " << *seqs_start << std::endl;


        bliss::index::kmer::KmerParser<KmerType > kmer_parser(fdata.valid_range_bytes);
        std::vector<KmerType> result;

        bool localcomp = true, comp = true;

        std::vector<unsigned char> gold;

        ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);

      //== loop over the reads
      for (; seqs_start != seqs_end; ++seqs_start)
      {
      if (fdata.valid_range_bytes.start <= (*seqs_start).id.get_pos()) {
        ++(this->seqCount);

              result.clear();

//              std::cout << "rank " << comm.rank() << " sequence " << *seqs_start << std::endl;

              // generate kmers and save them
              kmer_parser(*seqs_start, emplace_iter);

              gold.clear();
              gold.resize(seqs_start->seq_size(), 0);
              this->readFilePOSIX(this->fileName, seqs_start->seq_global_offset(), seqs_start->seq_size(), gold.data());
  //            std::cout << "rank " << comm.rank() << " entire gold: ";
  //            for (size_t j = 0; j < seqs_start->seq_size(); ++j) {
  //                                std::cout << gold[j];
  //                              }
  //            std::cout << std::endl;

              // compare the results.
              for (size_t i = 0; i < result.size(); ++i) {

                localcomp = equal(bliss::utils::KmerUtils::toASCIIString(result[i]).begin(), &(gold[i]), kmer_size);

                if (!localcomp) {

                  std::cout << "rank " << comm.rank() << " kmer: " << bliss::utils::KmerUtils::toASCIIString(result[i]) << std::endl;
                  std::cout << "rank " << comm.rank() << " gold: ";
                  for (size_t j = 0; j < kmer_size; ++j) {
                    std::cout << gold[i+j];
                  }
                  std::cout << std::endl;
                }

                comp &= localcomp;
              }  // end kmers for

              this->kmerCount += result.size();

              ASSERT_TRUE(comp);

      }
    //            std::cout << *seqs_start << ", ";
    //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
      //      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
      }
  //    int rank;
  //    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //    std::cout << "rank " << rank << " elem count " << this->seqCount << std::endl;

  //    if (this->seqCount == 0) {
  //      std::cout << " mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
  //    }

  }

};


TEST_P(NoNFASTQParseTest, parse_mmap)
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

TEST_P(NoNFASTQParseTest, parse_stdio)
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

TEST_P(NoNFASTQParseTest, parse_posix)
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
TEST_P(NoNFASTQParseTest, parse_mmap_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(NoNFASTQParseTest, parse_stdio_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(NoNFASTQParseTest, parse_posix_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(NoNFASTQParseTest, parse_mpiio_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::mpiio_file<bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

    this->parse_mpi(fobj, 0UL, comm);

    comm.barrier();
}
#endif



// using 21-mer, first entry is number of records, second entry is file size.
INSTANTIATE_TEST_CASE_P(Bliss, NoNFASTQParseTest, ::testing::Values(
    TestFileInfo(243,    434,     27580, std::string("/test/data/natural.fastq")),
    TestFileInfo(245,    490,     29250, std::string("/test/data/natural.withN.fastq")),
    TestFileInfo(7,      182,     939, std::string("/test/data/test.debruijn.small.fastq")),
    TestFileInfo(1,      26,      134, std::string("/test/data/test.debruijn.tiny.fastq")),
//    TestFileInfo(254562, 6618612, 34111308, std::string("/test/data/test.fastq")),
    TestFileInfo(140,    3640,    18761, std::string("/test/data/test.medium.fastq")),
    TestFileInfo(7,      182,     939, std::string("/test/data/test.small.fastq")),
    TestFileInfo(1,      16552,   33194, std::string("/test/data/test.unitiq1.fastq")),
    TestFileInfo(1,      3527,    7144, std::string("/test/data/test.unitiq1.short2.fastq")),
    TestFileInfo(1,      7367,    14824, std::string("/test/data/test.unitiq1.short.fastq")),
    TestFileInfo(1,      14603,   29296, std::string("/test/data/test.unitiq2.fastq")),
    TestFileInfo(2,      31155,   62490, std::string("/test/data/test.unitiqs.fastq"))
));


class SplitFASTQParseTest : public KmerReaderTest
{
  // can't compare the sequences directly since offsets may be different.
protected:
  template <typename ITER>
  using ParserType = bliss::io::FASTQParser<ITER>;

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

      ASSERT_TRUE(fdata.data[0] == '@');

      this->seqCount = 0;
      this->kmerCount = 0;

      // from FileLoader type, get the block iter type and range type
      using BlockIterType = typename ::bliss::io::file_data::const_iterator;
//    using OrigSeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;
//    using SeqIterType = ::bliss::iterator::filter_iterator<::bliss::io::SequenceNPredicate, OrigSeqIterType>;
      using SeqIterType = ::bliss::io::NSplitSequencesIterator<BlockIterType, bliss::io::FASTQParser >;

      ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

      size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes);

        //==  and wrap the chunk inside an iterator that emits Reads.
//        SeqIterType seqs_start(OrigSeqIterType(l1parser, fdata.begin(),
//            fdata.in_mem_end(), offset), OrigSeqIterType(fdata.in_mem_end()));
//        SeqIterType seqs_end(OrigSeqIterType(fdata.in_mem_end()));
        SeqIterType seqs_start(l1parser, fdata.begin(), fdata.in_mem_end(), offset);
        SeqIterType seqs_end(fdata.in_mem_end());

          bliss::index::kmer::KmerParser<KmerType > kmer_parser(fdata.valid_range_bytes);
          std::vector<KmerType> result;

          bool localcomp = true, comp = true;

          std::vector<unsigned char> gold;

          ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);


        //== loop over the reads
        for (; seqs_start != seqs_end; ++seqs_start)
        {
          auto ss = *seqs_start;
        if (fdata.valid_range_bytes.start <= ss.id.get_pos()) {
          ++(this->seqCount);

                result.clear();

                // generate kmers and save them
                kmer_parser(ss, emplace_iter);

                gold.clear();
                gold.resize(ss.seq_size(), 0);
                this->readFilePOSIX(this->fileName, ss.seq_global_offset(), ss.seq_size(), gold.data());

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
        }

  }


  template <typename file_loader>
  void parse_mpi(file_loader & fobj, size_t const & overlap, mxx::comm const & comm) {

      // get this->fileName

    ::bliss::io::file_data fdata = fobj.read_file();

//    std::cout << "rank " << comm.rank() << fobj.get_class_name() << std::endl;

//    comm.barrier();
//        if (comm.rank() == 1) {
//        std::cout << "rank " << comm.rank() << "   parent " << fdata.parent_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   inmem " << fdata.in_mem_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   search " << fdata.valid_range_bytes << std::endl;
//
//       std::cout << "rank " << comm.rank() << " ";
//       for (auto ii = fdata.data.begin(); ii != (fdata.data.begin() + (fdata.in_mem_range_bytes.end - fdata.in_mem_range_bytes.start)); ++ii ) {
//         std::cout << static_cast<unsigned char>(*ii);
//       }
//       std::cout << std::endl;
//        }
//        comm.barrier();

//    std::cout << " rank " << comm.rank() << " in mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

    ASSERT_TRUE(fobj.size() > 0);

    // parent range should match.
    ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
    ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

    // get the ranges to makes sure they are in memory.
    ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
    ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

    ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
    ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);

    if (fdata.valid_range_bytes.size() > 0) {
      if (fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] != '@') {
        std::cout << "ERROR: first char should be @, but is " << static_cast<unsigned char>(fdata.data[0]) << std::endl;
      }
      ASSERT_TRUE(fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] == '@');
    }

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
//      if (comm.rank() == 0) std::cout << "end "  << (i-1) << ": " << ends[i-1] << " begin next " << i << ": " << begins[i] << " overlap = " << overlap << std::endl;
      if (begins[i] == ends[i]) continue;

      ASSERT_TRUE(ends[j] == (begins[i]) );  // valid ranges do not include overlap.
      j = i;
    }
  //    std::cout << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

      this->seqCount = 0;
      this->kmerCount = 0;

    // from FileLoader type, get the block iter type and range type
    using BlockIterType = typename ::bliss::io::file_data::const_iterator;
//    using OrigSeqIterType = ::bliss::io::SequencesIterator<BlockIterType, bliss::io::FASTQParser >;
//    using SeqIterType = ::bliss::iterator::filter_iterator<::bliss::io::SequenceNPredicate, OrigSeqIterType>;
    using SeqIterType = ::bliss::io::NSplitSequencesIterator<BlockIterType, bliss::io::FASTQParser >;

    ParserType<typename ::bliss::io::file_data::const_iterator> l1parser;

    size_t offset = l1parser.init_parser(fdata.in_mem_cbegin(), fdata.parent_range_bytes, fdata.in_mem_range_bytes, fdata.valid_range_bytes, comm);

      //==  and wrap the chunk inside an iterator that emits Reads.
//      SeqIterType seqs_start(OrigSeqIterType(l1parser, fdata.begin(),
//          fdata.in_mem_end(), offset), OrigSeqIterType(fdata.in_mem_end()));
//      SeqIterType seqs_end(OrigSeqIterType(fdata.in_mem_end()));
    SeqIterType seqs_start(l1parser, fdata.begin(), fdata.in_mem_end(), offset);
    SeqIterType seqs_end(fdata.in_mem_end());


//          std::cout << "rank " << comm.rank() << " offset = " << offset << " sequence " << *seqs_start << std::endl;


        bliss::index::kmer::KmerParser<KmerType > kmer_parser(fdata.valid_range_bytes);
        std::vector<KmerType> result;

        bool localcomp = true, comp = true;

        std::vector<unsigned char> gold;

        ::fsc::back_emplace_iterator<std::vector<KmerType> > emplace_iter(result);

      //== loop over the reads
      for (; seqs_start != seqs_end; ++seqs_start)
      {
        auto ss = *seqs_start;
      if (fdata.valid_range_bytes.start <= ss.id.get_pos()) {
        ++(this->seqCount);

              result.clear();

//              std::cout << "rank " << comm.rank() << " sequence " << *seqs_start << std::endl;

              // generate kmers and save them
              kmer_parser(ss, emplace_iter);

              gold.clear();
              gold.resize(ss.seq_size(), 0);
              this->readFilePOSIX(this->fileName, ss.seq_global_offset(), ss.seq_size(), gold.data());
  //            std::cout << "rank " << comm.rank() << " entire gold: ";
  //            for (size_t j = 0; j < seqs_start->seq_size(); ++j) {
  //                                std::cout << gold[j];
  //                              }
  //            std::cout << std::endl;

              // compare the results.
              for (size_t i = 0; i < result.size(); ++i) {

                localcomp = equal(bliss::utils::KmerUtils::toASCIIString(result[i]).begin(), &(gold[i]), kmer_size);

                if (!localcomp) {

                  std::cout << "rank " << comm.rank() << " kmer: " << bliss::utils::KmerUtils::toASCIIString(result[i]) << std::endl;
                  std::cout << "rank " << comm.rank() << " gold: ";
                  for (size_t j = 0; j < kmer_size; ++j) {
                    std::cout << gold[i+j];
                  }
                  std::cout << std::endl;
                }

                comp &= localcomp;
              }  // end kmers for

              this->kmerCount += result.size();

              ASSERT_TRUE(comp);

      }
    //            std::cout << *seqs_start << ", ";
    //            std::cout << std::distance((*seqs_start).seq_begin, (*seqs_start).seq_end) << std::endl;
      //      std::cout << std::distance(l1.begin(), (*seqs_start).qual_begin) << ", " << std::distance(l1.begin(), (*seqs_start).qual_end) << std::endl;
      }
  //    int rank;
  //    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //    std::cout << "rank " << rank << " elem count " << this->seqCount << std::endl;

  //    if (this->seqCount == 0) {
  //      std::cout << " mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;
  //    }

  }



  template <typename file_loader>
  void parse_seq_file_helper(file_loader & fobj ) {
      bliss::io::file_data fdata = fobj.read_file();

      ASSERT_TRUE(fobj.size() > 0);

      ASSERT_EQ(fdata.getRange().size(), fobj.size());
      ASSERT_EQ(fdata.getRange().end, fobj.size());

      ASSERT_EQ(fdata.in_mem_range_bytes.size(), fobj.size());
      ASSERT_EQ(fdata.in_mem_range_bytes.end, fobj.size());

      ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
      ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

      ASSERT_TRUE(fdata.data[0] == '@');


      std::vector<KmerType> result;

      std::tie(this->seqCount, this->kmerCount) = bliss::io::KmerFileHelper::parse_file_data<bliss::index::kmer::KmerParser<KmerType >,
        bliss::io::FASTQParser, ::bliss::io::NSplitSequencesIterator>(fdata, result);

  }


  template <typename file_loader>
  void parse_mpi_file_helper(file_loader & fobj, size_t const & overlap, mxx::comm const & comm) {

      // get this->fileName

    ::bliss::io::file_data fdata = fobj.read_file();

//    std::cout << "rank " << comm.rank() << fobj.get_class_name() << std::endl;

//    comm.barrier();
//        if (comm.rank() == 1) {
//        std::cout << "rank " << comm.rank() << "   parent " << fdata.parent_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   inmem " << fdata.in_mem_range_bytes << std::endl;
//        std::cout << "rank " << comm.rank() << "   search " << fdata.valid_range_bytes << std::endl;
//
//       std::cout << "rank " << comm.rank() << " ";
//       for (auto ii = fdata.data.begin(); ii != (fdata.data.begin() + (fdata.in_mem_range_bytes.end - fdata.in_mem_range_bytes.start)); ++ii ) {
//         std::cout << static_cast<unsigned char>(*ii);
//       }
//       std::cout << std::endl;
//        }
//        comm.barrier();

//    std::cout << " rank " << comm.rank() << " in mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

    ASSERT_TRUE(fobj.size() > 0);

    // parent range should match.
    ASSERT_EQ(fdata.parent_range_bytes.size(), fobj.size());
    ASSERT_EQ(fdata.parent_range_bytes.end, fobj.size());

    // get the ranges to makes sure they are in memory.
    ASSERT_TRUE(fdata.in_mem_range_bytes.start >= fdata.parent_range_bytes.start);
    ASSERT_TRUE(fdata.in_mem_range_bytes.end <= fdata.parent_range_bytes.end);

    ASSERT_TRUE(fdata.in_mem_range_bytes.start <= fdata.valid_range_bytes.start);
    ASSERT_TRUE(fdata.in_mem_range_bytes.end >= fdata.valid_range_bytes.end);

    if (fdata.valid_range_bytes.size() > 0) {
      if (fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] != '@') {
        std::cout << "ERROR: first char should be @, but is " << static_cast<unsigned char>(fdata.data[0]) << std::endl;
      }
      ASSERT_TRUE(fdata.data[fdata.valid_range_bytes.start - fdata.in_mem_range_bytes.start] == '@');
    }

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
//      if (comm.rank() == 0) std::cout << "end "  << (i-1) << ": " << ends[i-1] << " begin next " << i << ": " << begins[i] << " overlap = " << overlap << std::endl;
      if (begins[i] == ends[i]) continue;

      ASSERT_TRUE(ends[j] == (begins[i]) );  // valid ranges do not include overlap.
      j = i;
    }
  //    std::cout << " orig mem " << fdata.in_mem_range_bytes << " valid " << fdata.valid_range_bytes << std::endl;

    std::vector<KmerType> result;

    std::tie(this->seqCount, this->kmerCount) = bliss::io::KmerFileHelper::parse_file_data<bliss::index::kmer::KmerParser<KmerType >,
      bliss::io::FASTQParser, ::bliss::io::NSplitSequencesIterator>(fdata, result, comm);

  }


};


TEST_P(SplitFASTQParseTest, parse_mmap)
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

TEST_P(SplitFASTQParseTest, parse_stdio)
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

TEST_P(SplitFASTQParseTest, parse_posix)
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
TEST_P(SplitFASTQParseTest, parse_mmap_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(SplitFASTQParseTest, parse_stdio_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(SplitFASTQParseTest, parse_posix_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(SplitFASTQParseTest, parse_mpiio_mpi)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::mpiio_file<bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

    this->parse_mpi(fobj, 0UL, comm);

    comm.barrier();
}
#endif

TEST_P(SplitFASTQParseTest, parse_mmap_helper)
{
#ifdef USE_MPI
    ::mxx::comm comm;
  if (comm.rank() == 0) {
#endif

    bliss::io::mmap_file fobj(this->fileName);

    this->parse_seq_file_helper(fobj);

#ifdef USE_MPI
  }
#endif

}

TEST_P(SplitFASTQParseTest, parse_stdio_helper)
{
#ifdef USE_MPI
    ::mxx::comm comm;
  if (comm.rank() == 0) {
#endif

    bliss::io::stdio_file fobj(this->fileName);

    this->parse_seq_file_helper(fobj);

#ifdef USE_MPI
  }
#endif

}

TEST_P(SplitFASTQParseTest, parse_posix_helper)
{
#ifdef USE_MPI
    ::mxx::comm comm;
  if (comm.rank() == 0) {
#endif

    bliss::io::posix_file fobj(this->fileName);

    this->parse_seq_file_helper(fobj);

#ifdef USE_MPI
  }
#endif

}

#ifdef USE_MPI
TEST_P(SplitFASTQParseTest, parse_mmap_mpi_helper)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::mmap_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi_file_helper(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(SplitFASTQParseTest, parse_stdio_mpi_helper)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::stdio_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi_file_helper(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(SplitFASTQParseTest, parse_posix_mpi_helper)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::partitioned_file<::bliss::io::posix_file, bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

  this->parse_mpi_file_helper(fobj, 0UL, comm);

  comm.barrier();
}

TEST_P(SplitFASTQParseTest, parse_mpiio_mpi_helper)
{
    ::mxx::comm comm;

  ::bliss::io::parallel::mpiio_file<bliss::io::FASTQParser> fobj(this->fileName, 0UL, comm);

    this->parse_mpi_file_helper(fobj, 0UL, comm);

    comm.barrier();
}
#endif

// using 21-mer, first entry is number of records, second entry is file size.
INSTANTIATE_TEST_CASE_P(Bliss, SplitFASTQParseTest, ::testing::Values(
    TestFileInfo(243,    434,     27580, std::string("/test/data/natural.fastq")),
    TestFileInfo(256,    490,     29250, std::string("/test/data/natural.withN.fastq")),
    TestFileInfo(7,      182,     939, std::string("/test/data/test.debruijn.small.fastq")),
    TestFileInfo(1,      26,      134, std::string("/test/data/test.debruijn.tiny.fastq")),
//    TestFileInfo(254562, 6618612, 34111308, std::string("/test/data/test.fastq")),
    TestFileInfo(140,    3640,    18761, std::string("/test/data/test.medium.fastq")),
    TestFileInfo(7,      182,     939, std::string("/test/data/test.small.fastq")),
    TestFileInfo(1,      16552,   33194, std::string("/test/data/test.unitiq1.fastq")),
    TestFileInfo(1,      3527,    7144, std::string("/test/data/test.unitiq1.short2.fastq")),
    TestFileInfo(1,      7367,    14824, std::string("/test/data/test.unitiq1.short.fastq")),
    TestFileInfo(1,      14603,   29296, std::string("/test/data/test.unitiq2.fastq")),
    TestFileInfo(2,      31155,   62490, std::string("/test/data/test.unitiqs.fastq"))
));


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



