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

/*
 * de_bruijn_node_trait.hpp
 *
 *  Created on: Jul 31, 2015
 *      Author: yongchao
 */

#ifndef DE_BRUIJN_NODE_TRAIT_HPP_
#define DE_BRUIJN_NODE_TRAIT_HPP_

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

#include "utils/logging.h"
#include "common/alphabets.hpp"

namespace bliss
{
	namespace de_bruijn
	{
		namespace node
		{
      static constexpr unsigned char SENSE = 0;
      static constexpr unsigned char ANTI_SENSE = 1;

      // utilities operating on nodes
      template <typename Kmer, typename EdgeType>
      class node_utils {
        public:

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_out_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<Kmer> & neighbors) {
            neighbors.clear();

            for (int i = 0; i < 4; ++i) {
              if (edge.get_edge_frequency(i) > 0) {
                neighbors.emplace_back(kmer);
                neighbors.back().nextFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_in_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<Kmer> & neighbors) {
            neighbors.clear();

            for (int i = 0; i < 4; ++i) {
              if (edge.get_edge_frequency(i + 4) > 0) {
                neighbors.emplace_back(kmer);
                neighbors.back().nextReverseFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_out_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<std::pair<Kmer, typename EdgeType::CountType> > & neighbors) {
            neighbors.clear();

            typename EdgeType::CountType count;
            for (int i = 0; i < 4; ++i) {
              count = edge.get_edge_frequency(i);
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextFromChar(i);
              }
            }
          }

          // construct a new kmer from a known edge, if that edge's count is non-zero
          static void get_in_neighbors(Kmer const & kmer, EdgeType const & edge, std::vector<std::pair<Kmer, typename EdgeType::CountType> > & neighbors) {
            neighbors.clear();

            typename EdgeType::CountType count;

            for (int i = 0; i < 4; ++i) {
              count = edge.get_edge_frequency(i + 4);
              if (count > 0) {
                neighbors.emplace_back(kmer, count);
                neighbors.back().first.nextReverseFromChar(i);
              }
            }
          }


      };


      // utilities operating on input edge type
      class input_edge_utils {
        public:
          template <typename Alphabet>
          static uint8_t reverse_complement_edges(uint8_t const & exts) {
            return (Alphabet::TO_COMPLEMENT[exts & 0xF] << 4) | Alphabet::TO_COMPLEMENT[exts >> 4];
          }

          template <typename Alphabet>
          static uint16_t reverse_complement_edges(uint16_t const & exts) {
            return (bliss::common::DNA16::TO_ASCII[bliss::common::DNA16::TO_COMPLEMENT[bliss::common::DNA16::FROM_ASCII[exts & 0xFF]]] << 8) |
                bliss::common::DNA16::TO_ASCII[bliss::common::DNA16::TO_COMPLEMENT[bliss::common::DNA16::FROM_ASCII[exts >> 8]]];
          }
      };

			/**
			 * @brief kmer metadata holding the incoming and outgoing edge counts in a de bruijn graph.
			 * @details lower 4 integers are out edge counts, upper 4 integers are in edge counts.
			 *        the out and in edges are relative to the kmer's orientation.
			 *
			 */
			template<typename ALPHA, typename COUNT = uint32_t>
			class edge_counts {

			  protected:
			    // 32 bit and 64 bit do not clamp.
			    template <typename C = COUNT, typename std::enable_if<sizeof(C) < 4, int>::type = 0 >
			    void clamped_add_1(C &target) {
			      if (target < std::numeric_limits<C>::max()) {
			        ++target;
			      }
			    }

          template <typename C = COUNT, typename std::enable_if<sizeof(C) >= 4, int>::type = 0  >
          void clamped_add_1(C &target) {
            ++target;
          }
          // 32 bit and 64 bit do not clamp.
          template <typename C = COUNT, typename std::enable_if<sizeof(C) < 4, int>::type = 0 >
          void clamped_add_1(C &target, uint8_t flag) {
            if ((target < std::numeric_limits<C>::max()) && (flag > 0)) {
              ++target;
            }
          }
          template <typename C = COUNT, typename std::enable_if<sizeof(C) >= 4, int>::type = 0  >
          void clamped_add_1(C &target, uint8_t flag) {
            if (flag > 0) ++target;
          }


			  public:

          using Alphabet = ALPHA;
          using CountType = COUNT;

	        friend std::ostream& operator<<(std::ostream& ost, const edge_counts<ALPHA, COUNT> & node)
	        {
	          // friend keyword signals that this overrides an externally declared function
	          ost << " dBGr node: counts self = " << node.counts[node.counts.size() - 1] << " in = [";
	          for (int i = 4; i < 8; ++i) ost << node.counts[i] << ",";
	          ost << "], out = [";
	          for (int i = 0; i < 4; ++i) ost << node.counts[i] << ",";
	          ost << "]";
	          return ost;
	        }


			    /// array of counts.  format:  [out A C G T; in A C G T; kmer count], ordered for the canonical strand, not necessarily same as for the input kmer..
			    std::array<COUNT, 9> counts;

				/*constructor*/
				edge_counts() : counts({{0, 0, 0, 0, 0, 0, 0, 0, 0}}) {};

				/*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
				~edge_counts() {}

				/**
				 * @brief update the current count from the input left and right chars (encoded in exts).
				 * @param relative_strand
				 * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  the exts should be consistent
				 *                        with the associated kmer's orientation.
				 */
				void update(uint8_t exts)
				{
				  clamped_add_1(counts[8]);   // increment self count.

//        // REQUIRE THAT EXTS has same orientation as associated kmer.
//				  // shuffle if antisense
//				  if (relative_strand == ANTI_SENSE) {  // swap upper and lower 4 bits
//				    exts = (exts << 4) | (exts >> 4);
//				  }

				  if (std::is_same<ALPHA, bliss::common::DNA>::value ||
				      std::is_same<ALPHA, bliss::common::DNA5>::value ||
				      std::is_same<ALPHA, bliss::common::RNA>::value ||
				      std::is_same<ALPHA, bliss::common::RNA5>::value) {
				    // value encodes only 1 possible character.  assume values are 0 1 2 3 for ACGT

				    if (std::is_same<ALPHA, bliss::common::DNA5>::value || std::is_same<ALPHA, bliss::common::RNA5>::value) {
				    	if ((exts & 0x0C) == 0) clamped_add_1(counts[exts & 0x3]);  // right edge (out).  increment if it's not marked as unknown.
				    	if ((exts & 0xC0) == 0) clamped_add_1(counts[((exts >> 4) & 0x3) + 4]);  // left edge (in).  increment if it's not marked as unknown.
				    } else {
				      // no clipping necessary.
				      clamped_add_1(counts[exts & 0x3]);  // right edge (out)
				      clamped_add_1(counts[((exts >> 4) & 0x3) + 4]);  // left edge (in)
				    }

				  } else {
				    // value encodes 0 or more possible characters.  bit position are ACGT from low to high

				    // if not DNA 16, convert to DNA16
				    if (!std::is_same<ALPHA, bliss::common::DNA16>::value) {
				      exts = (bliss::common::DNA16::FROM_ASCII[ALPHA::TO_ASCII[exts >> 4]] << 4) | bliss::common::DNA16::FROM_ASCII[ALPHA::TO_ASCII[exts & 0xF]];
				    }

            // now increment the counts.  Follows 4 bit bit ordering of DNA16, i.e. ACGT from lowest to highest.
            for (int i = 0; i < 8; ++i) {
              clamped_add_1(counts[i], (exts >> i) & 1);
            }
				  }

				}

				/**
				 * @brief update the current count from the input left and right chars (not encoded in exts).
				 *
				 * @param relative_strand
				 * @param exts              2 byte chars, in [out, in] lower byte being out.  order for the input kmer (not necessarily same as for the canonical)
				 */
				void update(uint16_t exts)
				{
          clamped_add_1(counts[8]);   // increment self count.

          uint8_t temp = (bliss::common::DNA16::FROM_ASCII[exts >> 8] << 4) | bliss::common::DNA16::FROM_ASCII[exts & 0xFF];

          // now increment the counts.  Follows 4 bit bit ordering of DNA16, i.e. ACGT from lowest to highest.
          for (int i = 0; i < 8; ++i) {
            clamped_add_1(counts[i], (temp >> i) & 1);
          }
				}

				COUNT get_edge_frequency(uint8_t idx) const {
				  if (idx >= 8) return 0;

				  return counts[idx];
				}

			};


      /*node trait class*/
      template<typename ALPHA>
      class edge_exists {
        public:

          using Alphabet = ALPHA;
          using CountType = uint8_t;

          friend std::ostream& operator<<(std::ostream& ost, const edge_exists<ALPHA> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (int i = 4; i < 8; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (int i = 0; i < 4; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }


          /// array of flags.  bit set to 1 if edge exists.  order from low to high bit:  Out A C G T; In A C G T. DNA 16 encoding.
          uint8_t counts;

        /*constructor*/
        edge_exists() : counts(0) {};

        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
        ~edge_exists() {}

        /**
         *
         * @param relative_strand
         * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
         */
        void update(uint8_t exts)
        {
//            printf("ext value before:  %x\n", exts);
          // convert to DNA16 first
          if (!std::is_same<ALPHA, bliss::common::DNA16>::value) {
            exts = (bliss::common::DNA16::FROM_ASCII[ALPHA::TO_ASCII[exts >> 4]] << 4) | bliss::common::DNA16::FROM_ASCII[ALPHA::TO_ASCII[exts & 0x0F]];
          }
//          printf("ext value converted:  %x\n", exts);
          counts |= exts;
        }

        /**
         *
         *
         * @param relative_strand
         * @param exts              2 byte chars, in [out, in] lower byte being out.  order for the input kmer (not necessarily same as for the canonical)
         */
        void update(uint16_t exts)
        {

          // construct a 2x4bit char.  no reordering.
          uint8_t temp = (bliss::common::DNA16::FROM_ASCII[exts >> 8] << 4) | bliss::common::DNA16::FROM_ASCII[exts & 0xFF];

          // update counts.
          counts |= temp;
        }

        uint8_t get_edge_frequency(uint8_t idx) const {
          if (idx >= 8) return 0;

          return (counts >> idx) & 0x1;
        }


      };



//      /*define the strand*/
//      static constexpr unsigned char SENSE = 0;
//      static constexpr unsigned char ANTI_SENSE = 1;
//
//      /*node trait class*/
//      template<typename Alphabet, typename IntType = int32_t>
//      class node_trait{
//        public:
//          static constexpr size_t edges_size = (Alphabet::SIZE + 7) >> 3 << 1;
//          static constexpr size_t counts_size = Alphabet::SIZE * 2;
//
//      protected:
//        /*Each strand occupies full bytes. The bytes with lower indices
//         * are for the sense trand and the rest for the antisense strand*/
//        uint8_t edges[edges_size];
//
//        /*coverage of the edge*/
//        IntType edge_cov[counts_size];
//
//        /*multiplicity of the k-mer node*/
//        IntType node_multiplicity;
//
//      public:
//
//        /*constructor*/
//        node_trait(){
//          /*clear each edge*/
//          for(int32_t i = 0; i < edges_size; ++i){
//            edges[i] = 0;
//          }
//          /*clear the coverage of each edge*/
//          for(int32_t i = 0; i < counts_size; ++i){
//            edge_cov[i] = 0;
//          }
//
//          /*initialize node multiplicity*/
//          node_multiplicity = 0;
//        }
//        /*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
//        ~node_trait()
//        {
//          /*do nothing*/
//        }
//
//        /*update node*/
//        void update(int32_t relative_strand, uint8_t exts)
//        {
//          /*an edge uses another encoding*/
//          uint8_t ch;
//          int32_t left_ext = -1, right_ext = -1;
//
//          /*decode the edge*/
//          if((ch = (exts >> 4) & 0x0f)){
//            /*from DNA16 to ASCII*/
//            ch = bliss::common::DNA16::TO_ASCII[ch];
//
//            /*from ASCII to DNA*/
//            left_ext = Alphabet::FROM_ASCII[ch];
//          }
//          if((ch = exts & 0x0f)){
//            /*from DNA16 to ASCII*/
//            ch = bliss::common::DNA16::TO_ASCII[ch];
//
//            /*from ASCII to DNA*/
//            right_ext = Alphabet::FROM_ASCII[ch];
//          }
//
//          /*update node trait*/
//          if(relative_strand == SENSE){
//            /*the input k-mer is identical to the cannonical k-mer*/
//            add_edge(SENSE, right_ext);
//            add_edge(ANTI_SENSE, left_ext >= 0 ? Alphabet::TO_COMPLEMENT[left_ext] : -1);
//          }else{
//            /*the input k-mer is not identical to the cannomical k-mer*/
//            add_edge(SENSE, left_ext >= 0 ? Alphabet::TO_COMPLEMENT[left_ext] : -1);
//            add_edge(ANTI_SENSE, right_ext);
//          }
//        }
//
//        void update(int32_t relative_strand, uint16_t exts)
//        {
//          uint8_t ch;
//          int32_t left_ext = -1, right_ext = -1;
//
//          /*decode the edge*/
//          if((ch = (exts >> 8) & 0x0ff)){
//            left_ext = Alphabet::FROM_ASCII[ch];
//          }
//          if((ch = exts & 0x0ff)){
//            right_ext = Alphabet::FROM_ASCII[ch];
//          }
//
//          /*update node trait*/
//          if(relative_strand == SENSE){
//            /*the input k-mer is identical to the cannonical k-mer*/
//            add_edge(SENSE, right_ext);
//            add_edge(ANTI_SENSE, left_ext >= 0 ? Alphabet::TO_COMPLEMENT[left_ext] : -1);
//          }else{
//            /*the input k-mer is not identical to the cannomical k-mer*/
//            add_edge(SENSE, left_ext >= 0 ? Alphabet::TO_COMPLEMENT[left_ext] : -1);
//            add_edge(ANTI_SENSE, right_ext);
//          }
//        }
//
//        /*get the multiplicity of the k-mer node*/
//        inline IntType get_node_multiplicity(){
//          return node_multiplicity;
//        }
//
//        /**
//         * @brief get the number of edges (sum of counts) given the strand.  this is the "out" edges given the SENSE
//         */
//        IntType get_num_edges(int32_t strand){
//          /*simple strand test*/
//          assert(strand == SENSE || strand == ANTI_SENSE);
//
//          /*get the base address for the strand*/
//          IntType *ptr = edge_cov + strand * Alphabet::SIZE;
//
//          /*compute the sum*/
//          IntType sum = 0;
//          for(int i = 0; i < Alphabet::SIZE; ++i){
//            sum += ptr[i];
//          }
//          return sum;
//        }
//
//        /**
//         * @brief get the number of unique edges given the strand..  this is the "out" edges given the SENSE.
//         */
//        IntType get_num_unique_edges(int32_t strand){
//          /*simple strand test*/
//          assert(strand == SENSE || strand == ANTI_SENSE);
//
//          /*get the base address for the strand*/
//          IntType *ptr = edge_cov + strand * Alphabet::SIZE;
//
//          /*compute the sum*/
//          IntType sum = 0;
//          for(int i = 0; i < Alphabet::SIZE; ++i){
//            sum += (ptr[i] > 0 ? 1 : 0);
//          }
//          return sum;
//        }
//
//
//        /*add an edge given the strand and the base extension*/
//        void add_edge(int32_t strand, int32_t ext){
//          /*simple strand test*/
//          assert(strand == SENSE || strand == ANTI_SENSE);
//
//          /*simple base extension test*/
//          assert(ext >= 0 && ext < Alphabet::SIZE);
//
//
//          /*must be an effective edge*/
//          if(ext >= 0){
//            /*get the base address for the strand*/
//            uint8_t* ptr = edges + strand * ((Alphabet::SIZE + 7) >> 3);
//            IntType* ptr2 = edge_cov + strand * Alphabet::SIZE;
//
//            /*set the bit*/
//            ptr += ext >> 3;
//            *ptr |= 1 << (ext & 7);
//
//            /*increase edge coverage*/
//            ptr2[ext]++;
//          }
//          /*increase node multiplicity*/
//          ++node_multiplicity;
//        }
//
//        /*remove the edge*/
//        void remove_edge(int32_t strand, int32_t ext){
//          /*simple strand test*/
//          assert(strand == SENSE || strand == ANTI_SENSE);
//
//          /*simple base extension test*/
//          assert(ext >= 0 && ext < Alphabet::SIZE);
//
//
//          /*must be an effective edge*/
//           if(ext < 0){
//             cerr << "The extension base for the edge is invalid: " << ext << endl;
//             return;
//           }
//
//          /*get the base address for the strand*/
//          uint8_t* ptr = edges + strand * ((Alphabet::SIZE + 7) >> 3);
//          IntType* ptr2 = edge_cov + strand * Alphabet::SIZE;
//
//          /*simple strand test*/
//          assert(strand != SENSE && strand != ANTI_SENSE);
//
//          /*simple base extension test*/
//          assert(ext >= 0 && ext < Alphabet::SIZE);
//
//          /*clear the bit*/
//          ptr += ext >> 3;
//          *ptr &= ~(1 << (ext & 7));
//
//          /*clear coverage*/
//          ptr2[ext] = 0;
//        }
//
//
//        /*check if the edge exists*/
//        bool edge_exists(int32_t strand, int32_t ext){
//          /*simple strand test*/
//          assert(strand == SENSE || strand == ANTI_SENSE);
//
//          /*simple base extension test*/
//          assert(ext >= 0 && ext < Alphabet::SIZE);
//
//          /*must be an effective edge*/
//         if(ext < 0){
//           cerr << "The extension base for the edge is invalid: " << ext << endl;
//           return false;
//         }
//
//          /*get the base address for the strand*/
//          uint8_t* ptr = edges + strand * ((Alphabet::SIZE + 7) >> 3);
//
//
//
//          /*locate the byte address*/
//          ptr += ext >> 3;
//
//          return ((*ptr >> (ext & 7)) & 1) ? true : false;
//        }
//      };



		}/*namespace node*/
	}/*namespace de_bruijn*/
}/*namespace bliss*/




#endif /* DE_BRUIJN_NODE_TRAIT_HPP_ */
