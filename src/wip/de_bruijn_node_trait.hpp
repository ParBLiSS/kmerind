/*
 * de_bruijn_node_trait.hpp
 *
 *  Created on: Jul 31, 2015
 *      Author: yongchao
 */

#ifndef DE_BRUIJN_NODE_TRAIT_HPP_
#define DE_BRUIJN_NODE_TRAIT_HPP_
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
#include "common/base_types.hpp"

#include "utils/timer.hpp"

namespace bliss
{
	namespace de_bruijn
	{
		namespace node
		{
      static constexpr unsigned char SENSE = 0;
      static constexpr unsigned char ANTI_SENSE = 1;

			/*node trait class*/
			template<typename Alphabet, typename CountType = uint32_t>
			class edge_counts{
			  public:

	        friend std::ostream& operator<<(std::ostream& ost, const edge_counts<Alphabet, CountType> & node)
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
			    std::array<CountType, 9> counts;

				/*constructor*/
				edge_counts() : counts({{0, 0, 0, 0, 0, 0, 0, 0, 0}}) {};

				/*destructor.  not virtual, so that we don't have virtual lookup table pointer in the structure as well.*/
				~edge_counts() {}

				/**
				 *
				 * @param relative_strand
				 * @param exts            2 4bits in 1 uchar.  ordered as [out, in], lower bits being out.  ordered for the input kmer (not necessarily canonical)
				 */
				void update(uint8_t relative_strand, uint8_t exts)
				{
				  ++counts[8];   // increment self count.


				  // shuffle if antisense
				  if (relative_strand == ANTI_SENSE) {  // swap upper and lower 4 bits
				    exts = (exts << 4) | (exts >> 4);
				  }

          // now increment the counts.  Follows 4 bit bit ordering of DNA16, i.e. ACGT from lowest to highest.
				  uint8_t mask = 1;
				  for (int i = 0; i < 8; ++i) {
				    counts[i] += (exts >> i) & mask;
				  }

				}

				/**
				 *
				 *
				 * @param relative_strand
				 * @param exts              2 byte chars, in [out, in] lower byte being out.  order for the input kmer (not necessarily same as for the canonical)
				 */
				void update(uint8_t relative_strand, uint16_t exts)
				{

				  // construct a 2x4bit char.  no reordering.
				  uint8_t temp = bliss::common::DNA16::FROM_ASCII[exts & 0xFF] |
				        (bliss::common::DNA16::FROM_ASCII[exts >> 8] << 4);

				  // now increment the counts - delegate to other update() function.
				  update(relative_strand, temp);
				}

			};


      /*node trait class*/
      template<typename Alphabet>
      class edge_exists{
        public:

          friend std::ostream& operator<<(std::ostream& ost, const edge_exists<Alphabet> & node)
          {
            // friend keyword signals that this overrides an externally declared function
            ost << " dBGr node: in = [";
            for (int i = 4; i < 8; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "], out = [";
            for (int i = 0; i < 4; ++i) ost << (((node.counts >> i) & 0x1) == 1 ? 1 : 0) << ",";
            ost << "]";
            return ost;
          }


          /// array of counts.  format:  [out A C G T; in A C G T; kmer count], ordered for the canonical strand, not necessarily same as for the input kmer..
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
        void update(uint8_t relative_strand, uint8_t exts)
        {
          // shuffle if antisense
          if (relative_strand == ANTI_SENSE) {  // swap upper and lower 4 bits
            exts = (exts << 4) | (exts >> 4);
          }

          counts |= exts;
        }

        /**
         *
         *
         * @param relative_strand
         * @param exts              2 byte chars, in [out, in] lower byte being out.  order for the input kmer (not necessarily same as for the canonical)
         */
        void update(uint8_t relative_strand, uint16_t exts)
        {

          // construct a 2x4bit char.  no reordering.
          uint8_t temp = bliss::common::DNA16::FROM_ASCII[exts & 0xFF] |
                (bliss::common::DNA16::FROM_ASCII[exts >> 8] << 4);

          // now increment the counts - delegate to other update() function.
          update(relative_strand, temp);
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
