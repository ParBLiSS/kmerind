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
 * @file    bitgroup_ops.hpp
 * @ingroup
 * @author  tpan
 * @brief   collection of bit reverse routines for single word or a byte array (or array of words)
 * @details support reverse by bit groups, i.e. within bit group, the bit order does not change.
 *          Designed as a set of functors because this allows us to do partial template specialization (BIT GROUP SIZE and WORD_TYPE, whereas functions need to be fully specialized.
 *
 *          Additionally, the functions are limited in the following ways:
 *          1. input needs to be a word (32 or 64bit for SWAR, 128 bit for SSSE3, and 256 bit for AVX2).
 *          1.1.  if less than a word's worth, the underlying array is read past end, reversed, and result is memcpy'd while respecting boundaries.
 *          1.2.  if necessary to OR bytes together, the caller of the functors should do that.
 *          2. BIT_GROUP_SIZE should not be larger than word size.  (e.g. 8 2byte groups in 128 bits can be swapped using a single shuffle, excluding
 *              load and store.  should be faster than iterating linearly across the 16 bytes)
 *          3. bit offset, where necessary (e.g. when using non-power of 2 bitgroups), is limited to less than 8  (smaller than 1 byte)
 *          3.1.  if larger offset is needed, calling function should do the "shifting" of the output ptr to the floor of the shift in bytes.
 *          4. offset is a template parameter.  the calling function should properly dispatch to the right function instance.
 *          5. non-simd 64 bit shift probably sufficient - number of machine words to shift is known so compiler can unroll loop.
 *             SSSE3 and AVX2 do not have an advantage (6 instruction and 12 instructions, respectively.)
 *          6. use 64bit.  on 32 bit processors, let hardware emulate.
 *
 */
#ifndef SRC_UTILS_BITGROUP_OPS_HPP_
#define SRC_UTILS_BITGROUP_OPS_HPP_

#define BITS_INLINE inline

#include <type_traits>   // is_integral, etc
#include <limits>        // numeric_limits
#include <string>
#include <cstring>       // memcpy and memset
#include <sstream>  // stringstream
#include <cstdint>       // uint8_t

//#include <atomic>    // for ensuring memory ordering.

#include "utils/logging.h"
#include "bliss-config.hpp"
#if defined(__AVX2__) || defined(__SSSE3__)
#include <x86intrin.h>   // all intrinsics.  will be enabled based on compiler flag such as __SSSE3__ internally.
#endif

#if defined __GNUC__ && __GNUC__>=6
// disable __m128i and __m256i ignored attribute warning in gcc
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

// needed by clang else it does not know where to get the bswap function.  however, conflicts with farmhash.cc
//#include <byteswap.h>



// done:  3bit reverse - done for SWAR, SSSE3, AVX2
// done:  see effect of not shifting when working with byte arrays.  okay not to shift, if bit_offset is passed to reverse
// done:  vector reverse.
// TODO:  perf compare to old impl - slower.  cause:  branching, and sometimes non-inlining.
// TODO:  add other operations such as masking, shifting, loading data from mem and storing data to mem.  refactor array version of reverse to use load/op/store

// NOTE: alignas placement is accepted by GCC at end of array declaration, or between array name and size specification.  "name[size] alignas(bits)"  or "name alignas(bits) [size]"
//       clang accepts the second form only.
//        also, clang requires that variable definition to have the same alignment specification as in the declaration.


#if defined(HAVE_BUILTIN_BSWAP) || defined(__clang__) ||                \
  (defined(__GNUC__) && ((__GNUC__ == 4 && __GNUC_MINOR__ >= 8) ||      \
                         __GNUC__ >= 5))
#define _bliss_bswap_64(x)  __builtin_bswap64(x)
#define _bliss_bswap_32(x)  __builtin_bswap32(x)
#elif defined(__ICC)
#define _bliss_bswap_64(x)  _bswap64(x)
#define _bliss_bswap_32(x)  _bswap(x)
#endif

#define _bliss_bswap_16(x)  (((x >>  8) &  0xFF) | ((x &  0xFF) <<  8))

#define _bliss_abs(x)  ((x < 0) ? -x : x)

// TODO: reverse_groups vs reverse_in_group
//#define MEM_ORDER ::std::memory_order_seq_cst

namespace bliss {

  namespace utils {

    namespace bit_ops {


      static constexpr unsigned char BIT_REV_SEQ = 0;
      static constexpr unsigned char BIT_REV_SWAR = 1;   // SIMD Within A Register
      static constexpr unsigned char BIT_REV_SSSE3 = 2;
      static constexpr unsigned char BIT_REV_AVX2 = 4;

      // TODO: replace all unsigned char SIMD specifiers.
      // for now, this is used only for the generic operator version of reverse.
      struct BITREV_SEQ {
    	  using MachineWord = uint64_t;
    	  static constexpr unsigned char SIMDVal = BIT_REV_SEQ;
      };
      struct BITREV_SWAR {
    	  using MachineWord = uint64_t;
          static constexpr unsigned char SIMDVal = BIT_REV_SWAR;
      };
      struct BITREV_SSSE3 {
#if defined(__SSSE3__)
    	  using MachineWord = __m128i;
          static constexpr unsigned char SIMDVal = BIT_REV_SSSE3;
#else
    	  using MachineWord = uint64_t;
          static constexpr unsigned char SIMDVal = BIT_REV_SWAR;
#endif
      };

      struct BITREV_AVX2 {
#if defined(__AVX2__)
    	  using MachineWord = __m256i;
          static constexpr unsigned char SIMDVal = BIT_REV_AVX2;
#elif defined(__SSSE3__)
    	  using MachineWord = __m128i;
          static constexpr unsigned char SIMDVal = BIT_REV_SSSE3;
#else
    	  using MachineWord = uint64_t;
          static constexpr unsigned char SIMDVal = BIT_REV_SWAR;
#endif
      };

      /// automatically choose the most appropriate MachineWord and SIMD type based
      /// on supported simd capability and number of bytes to be reversed.
      template <size_t BYTES, typename MAX_SIMD = BITREV_AVX2>
      struct BITREV_AUTO_AGGRESSIVE {
    	  using MachineWord =
#if defined(__AVX2__)
    			  typename ::std::conditional<((BYTES > 16) && (MAX_SIMD::SIMDVal == BIT_REV_AVX2)), __m256i,
#endif
#if defined(__SSSE3__)
    			    typename ::std::conditional<((BYTES > 8) && (MAX_SIMD::SIMDVal >= BIT_REV_SSSE3)),  __m128i,
#endif
    			    typename ::std::conditional<(BYTES > 4),  uint64_t,
    			        typename ::std::conditional<(BYTES > 2),  uint32_t,
    			          typename ::std::conditional<(BYTES > 1),  uint16_t,
    			          	uint8_t
    			          >::type
    			        >::type
  			          >::type
#if defined(__SSSE3__)
  			        >::type
#endif
#if defined(__AVX2__)
		          >::type
#endif
    	  ;

    	  static constexpr unsigned char SIMDVal =
    			  (BYTES > 16) ?
    					MAX_SIMD::SIMDVal :
    					(BYTES > 8) ?
    					  ((MAX_SIMD::SIMDVal == BIT_REV_AVX2) ? BIT_REV_SSSE3 : MAX_SIMD::SIMDVal) :
    					  BIT_REV_SWAR;
      };

      /// automatically choose the most appropriate MachineWord and SIMD type based
      /// on supported simd capability and number of bytes to be reversed.
      template <size_t BYTES, typename MAX_SIMD = BITREV_AVX2>
      struct BITREV_AUTO_CONSERVATIVE {
    	  using MachineWord =
#if defined(__AVX2__)
    			  typename ::std::conditional<((BYTES >= 32) && (MAX_SIMD::SIMDVal == BIT_REV_AVX2)), __m256i,
#endif
#if defined(__SSSE3__)
    			    typename ::std::conditional<((BYTES >= 16) && (MAX_SIMD::SIMDVal >= BIT_REV_SSSE3)),  __m128i,
#endif
    			    typename ::std::conditional<(BYTES >= 8),  uint64_t,
    			        typename ::std::conditional<(BYTES >= 4),  uint32_t,
    			          typename ::std::conditional<(BYTES >= 2),  uint16_t,
    			          	uint8_t
    			          >::type
    			        >::type
  			          >::type
#if defined(__SSSE3__)
  			        >::type
#endif
#if defined(__AVX2__)
		          >::type
#endif
    	  ;

    	  static constexpr unsigned char SIMDVal =
    			  (BYTES >= 32) ?
    					MAX_SIMD::SIMDVal :
    					(BYTES >= 16) ?
    					  ((MAX_SIMD::SIMDVal == BIT_REV_AVX2) ? BIT_REV_SSSE3 : MAX_SIMD::SIMDVal) :
    					  BIT_REV_SWAR;
      };




      /// shift right by number of bits less than 8.
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE
      srli(WORD_TYPE const & val, uint16_t shift) {
        return val >> shift;
      }
      /// shift right by number of bits less than 8.
      template <uint16_t SHIFT, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(SHIFT == 0), WORD_TYPE>::type
      srli(WORD_TYPE const & val) {
        return val;
      }
      template <uint16_t SHIFT, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(SHIFT >= (sizeof(WORD_TYPE) << 3)), WORD_TYPE>::type
      srli(WORD_TYPE const & val) {
        return 0;
      }
      template <uint16_t SHIFT, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(SHIFT > 0) && (SHIFT < (sizeof(WORD_TYPE) << 3)), WORD_TYPE>::type
      srli(WORD_TYPE const & val) {
        return val >> SHIFT;
      }

      /// shift left by number of bits less than 8.
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE
      slli(WORD_TYPE const & val, uint16_t shift) {
        return val << shift;
      }
      /// shift left by number of bits less than 8.
      template <uint16_t SHIFT, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(SHIFT == 0), WORD_TYPE>::type
      slli(WORD_TYPE const & val) {
        return val;
      }
      template <uint16_t SHIFT, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(SHIFT >= (sizeof(WORD_TYPE) << 3)), WORD_TYPE>::type
      slli(WORD_TYPE const & val) {
        return 0;
      }
      template <uint16_t SHIFT, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(SHIFT > 0) && (SHIFT < (sizeof(WORD_TYPE) << 3)), WORD_TYPE>::type
      slli(WORD_TYPE const & val) {
        return val << SHIFT;
      }
// these cause tests to fail. explicit calls to slli and srli works, however.
//      template <int16_t SHIFT, typename WORD_TYPE>
//      BITS_INLINE typename ::std::enable_if<(SHIFT == 0), WORD_TYPE>::type
//      shift(WORD_TYPE const & val) {
//        return val;
//      }
//      template <int16_t SHIFT, typename WORD_TYPE>
//      BITS_INLINE typename ::std::enable_if<(SHIFT > 0), WORD_TYPE>::type
//      shift(WORD_TYPE const & val) {
//        return slli<SHIFT>(val);
//      }
//      template <int16_t SHIFT, typename WORD_TYPE>
//      BITS_INLINE typename ::std::enable_if<(SHIFT < 0), WORD_TYPE>::type
//      shift(WORD_TYPE const & val) {
//        return srli<(-SHIFT)>(val);
//      }

      template <typename DEST_WORD_TYPE, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<(sizeof(DEST_WORD_TYPE) <= 8), DEST_WORD_TYPE>::type
      loadu(WORD_TYPE const * u) {
        return *(reinterpret_cast<DEST_WORD_TYPE const *>(u));
      }

      template <typename SRC_WORD_TYPE, typename WORD_TYPE,
      typename = typename ::std::enable_if<(sizeof(SRC_WORD_TYPE) <= 8)>::type>
      BITS_INLINE void storeu(WORD_TYPE * u, SRC_WORD_TYPE const & val) {
    	  *(reinterpret_cast<SRC_WORD_TYPE *>(u)) = val;
      }

      template <size_t nbytes, uint8_t byte_offset = 0, typename DEST_WORD_TYPE, typename WORD_TYPE,
       	   typename ::std::enable_if<((nbytes + byte_offset) <= sizeof(DEST_WORD_TYPE)), int>::type = 1>
      BITS_INLINE void load_part(DEST_WORD_TYPE & x, WORD_TYPE const *u) {
    	  memcpy(reinterpret_cast<uint8_t*>(&x) + byte_offset, u, nbytes);
      }

      template <size_t nbytes, uint8_t byte_offset = 0, typename WORD_TYPE, typename SRC_WORD_TYPE,
      	  typename ::std::enable_if<((nbytes + byte_offset) <= sizeof(SRC_WORD_TYPE)), int>::type = 1>
      BITS_INLINE void store_part(WORD_TYPE *u, SRC_WORD_TYPE const & val) {
    	  // TODO: faster way of doing this.  right now, this translates to 6 memory instructions
    	  // "movq".  could potentially be done via 2.
    	  memcpy(u, reinterpret_cast<uint8_t const *>(&val) + byte_offset, nbytes);
      }



      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE bit_not(WORD_TYPE const & u) {
        return ~u;
      }

      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE bit_or(WORD_TYPE const & u, WORD_TYPE const & v) {
        return u | v;
      }
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE bit_and(WORD_TYPE const & u, WORD_TYPE const & v) {
        return u & v;
      }
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE bit_xor(WORD_TYPE const & u, WORD_TYPE const & v) {
        return u ^ v;
      }
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE equal(WORD_TYPE const & u, WORD_TYPE const & v) {
        return u == v;
      }
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE less(WORD_TYPE const & u, WORD_TYPE const & v) {
        return u < v;
      }

      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE zero() {
        return 0;
      }
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE bit_max() {
        return bit_not(zero<WORD_TYPE>());
      }


      /**       * @brief base bit reverse type. base template only
       * @tparam WORD_TYPE          type of input word
       * @tparam BIT_GROUP_SIZE     number of bits in a group to be reversed.  supports any value less than 8, and powers of 2.  tested 3 and powers of 2.  cannot exceed word size.
       * @tparam BIT_REV_SIMD_TYPE  type of algorithm to use based on available hardware
       * @tparam DUMMY              here so that we can define static constexpr array variables in the header file.
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char BIT_REV_SIMD_TYPE = BIT_REV_SEQ, bool POW2 = ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0)>
      struct bitgroup_ops {

          //  default imple is for SEQuential bit reverse.  this is defined for all bit_group_sizes.
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < (sizeof(uint64_t) << 3), "ERROR: BIT_GROUP_SIZE is greater than number of bits in uint64_t");
          static_assert(((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0) || (BIT_GROUP_SIZE < 8), "ERROR: BIT_GROUP_SIZE is has to be powers of 2 or less than 8");

          /* Linear (inefficient) reverse: */

          static constexpr uint64_t group_mask = ~(std::numeric_limits<uint64_t>::max() << BIT_GROUP_SIZE);

          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_SIMD_TYPE;

          /**
           * @brief reverse function to reverse bits for data types are are not power of 2 number of bytes
           * @details OR the current len size result with any existing output, so important to make sure the underlying out vector is appropriately initialized.
           *          len has to be a multiple of BIT_GROUP_SIZE if BIT_GROUP_SIZE is greater than 8 (i.e. multiple bytes).
           *          in other words, when a converted word is divided into BIT_GROUP_SIZE, there should not be padding zero's in the groups
           *
           * @param len       length of input in bytes
           * @param out       pointer into a byte array.  the array should be initialized to 0 before any call to this function.
           * @return number of bits in the last byte that are not swapped (i.e. last 2 bits for 3 bits op.)
           */
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if ((len << 3) < BIT_GROUP_SIZE) return bit_offset; //throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            // if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be less than 8.  for larger, shift input instead.");
            assert(bit_offset < 8);

            if((len % ((BIT_GROUP_SIZE + 7) >> 3)) != 0) {
              printf("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes. len %lu, BITS %u", len, BIT_GROUP_SIZE);
            }
            assert((len % ((BIT_GROUP_SIZE + 7) >> 3)) == 0);

            // SHIFTS are NOT NECESSARY.

            switch (len) {
              case 8:
                {
                  uint64_t v8 = *(reinterpret_cast<const uint64_t*>(in));
                  // v8 >>= bit_offset;
                  *(reinterpret_cast<uint64_t*>(out)) |= this->reverse(v8, bit_offset); // >> bit_offset;
                }
                break;
              case 7:
              case 6:
              case 5:
                {
                  uint64_t v7i =  *(reinterpret_cast<const uint64_t*>(in));
                  //memcpy(reinterpret_cast<uint8_t*>(&v7), in, len);
                  ///v7 >>= bit_offset;
                  uint64_t v7 = this->reverse(v7i, bit_offset); // >> bit_offset;
                  out[0] |= *(reinterpret_cast<uint8_t*>(&v7) + 8 - len);
                  memcpy(out + 1, reinterpret_cast<uint8_t*>(&v7) + (9 - len), len - 2);
                  out[len - 1] |= *(reinterpret_cast<uint8_t*>(&v7) + 7);

                }
                break;
              case 4:
                {
                  uint32_t v4 = *(reinterpret_cast<const uint32_t*>(in));
                  //v4 >>= bit_offset;
                  *(reinterpret_cast<uint32_t*>(out)) |= this->reverse(v4, bit_offset); // >> bit_offset;
                }
                break;
              case 3:
                {
                  uint32_t v3i = *(reinterpret_cast<const uint32_t*>(in));
                  //memcpy(reinterpret_cast<uint8_t*>(&v3), in, 3);
                  // v3 >>= bit_offset;
                  uint32_t v3 = this->reverse(v3i, bit_offset); // >> bit_offset;
                  out[0] |= *(reinterpret_cast<uint8_t*>(&v3) + 1);
                  out[1] =  *(reinterpret_cast<uint8_t*>(&v3) + 2);
                  out[2] |= *(reinterpret_cast<uint8_t*>(&v3) + 3);

                }
                break;
              case 2:
                {
                  uint16_t v2 = *(reinterpret_cast<const uint16_t*>(in));
                  //v2 >>= bit_offset;
                  *(reinterpret_cast<uint16_t*>(out)) |= this->reverse(v2, bit_offset); // >> bit_offset;
                }
                break;
              case 1:
                {
                  uint8_t v1 = *in;
                  // v1 >>= bit_offset;
                  *out |= this->reverse(v1, bit_offset); // >> bit_offset;
                }
                break;
              default:
                BL_ERRORF("ERROR:  length of array should be > 0 and <= 8.  For longer, break it up or use the vector version of bitgroup_ops.");
                break;
            }

            return ((len << 3) - bit_offset) % BIT_GROUP_SIZE;
          }

          /**
           * @brief  reverse the lower portion of WORD, starting from 0, up to the largest multiples of BIT_GROUP_SIZE that still fits.
           *
           * @param v    input word.  LSB aligned to start of a bit group.
           * @return     bit reversed word.  MSB aligned to end of the originally lowest bit group.
           */
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE, typename std::enable_if<(BITS >= (sizeof(WORD_TYPE) << 3)), int>::type = 1>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & v, uint16_t bit_offset) const {
//            static_assert(BIT_GROUP_SIZE == (sizeof(WORD_TYPE) << 3), "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");
            return v;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE,
        		  typename std::enable_if<(BITS < (sizeof(WORD_TYPE) << 3)), int>::type = 1>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & v, uint16_t bit_offset) const {

            // copy the source data
            WORD_TYPE w = 0;
            WORD_TYPE u = v;

            unsigned int rem = ((sizeof(WORD_TYPE) << 3) - bit_offset) % BIT_GROUP_SIZE;
            //WORD_TYPE rem_mask = static_cast<WORD_TYPE>(~(std::numeric_limits<WORD_TYPE>::max() >> rem));

            u >>= bit_offset;

            // now do the middle part
            for (size_t i = bit_offset; i < ((sizeof(WORD_TYPE) << 3) - rem); i += BIT_GROUP_SIZE) {
              w <<= BIT_GROUP_SIZE;
              w |= (u & group_mask);
              u >>= BIT_GROUP_SIZE;
            }

            // don't put back the remainder.  just shift.. we'll need to OR bytes together anyways.
            w <<= rem;
//            w |= (u & rem_mask);

            // finally merge back to v.
            return w;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE, uint16_t bit_offset = 0, typename WORD_TYPE>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & v) const {
            return reverse(::std::forward<WORD_TYPE const>(v), bit_offset);
          }

      };

      /// partial template specialization for SWAR based bit reverse.  uses BSWAP when appropriate, else use bit shift  this is defined only for bit_group_sizes that are powers of 2 up to sizeof WORD_TYPE)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SWAR, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < (sizeof(uint64_t) << 3), "ERROR: BIT_GROUP_SIZE is greater than number of bits in uint64_t");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");


          // masks.  if uint64_t indicates that system is 32 bit, then this should truncate.
          static constexpr uint64_t mask32 = 0x00000000FFFFFFFF;
          static constexpr uint64_t mask16 = 0x0000FFFF0000FFFF;
          static constexpr uint64_t  mask8 = 0x00FF00FF00FF00FF;
          static constexpr uint64_t  mask4 = 0x0F0F0F0F0F0F0F0F;
          static constexpr uint64_t  mask2 = 0x3333333333333333;
          static constexpr uint64_t  mask1 = 0x5555555555555555;

          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_SWAR;

          /// reverse function to reverse bits for data types are are not power of 2 number of bytes
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if ( (len << 3) < BIT_GROUP_SIZE ) return bit_offset; // throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be 0, since we are ORing the first and last bytes.
            //if (bit_offset > 0) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be 0 for SWAR with power of 2 BIT_GROUP_SIZE.");
            assert(bit_offset == 0 );

            //if ((len % ((BIT_GROUP_SIZE + 7) >> 3)) > 0)
            //  throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
            assert((len % ((BIT_GROUP_SIZE + 7) >> 3)) == 0);


            switch (len) {
              case 8:
                *(reinterpret_cast<uint64_t*>(out)) = this->reverse(*(reinterpret_cast<const uint64_t*>(in)));
                break;
              case 7:
              case 6:
              case 5:
                {
                  uint64_t v64i = *(reinterpret_cast<const uint64_t*>(in));
                  //memcpy(reinterpret_cast<uint8_t*>(&v64), in, len);
                  uint64_t v64 = this->reverse(v64i);
                  memcpy(out, reinterpret_cast<uint8_t*>(&v64) + (8-len), len);
                }
                break;
              case 4:
                *(reinterpret_cast<uint32_t*>(out)) = this->reverse(*(reinterpret_cast<const uint32_t*>(in)));
                break;
              case 3:
                {
                  uint32_t v32i = *(reinterpret_cast<const uint32_t*>(in));
                  //memcpy(reinterpret_cast<uint8_t*>(&v32), in, len);
                  uint32_t v32 = this->reverse(v32i);
                  memcpy(out, reinterpret_cast<uint8_t*>(&v32) + (4 - len), len);
                }
                break;
              case 2:
                *(reinterpret_cast<uint16_t*>(out)) = this->reverse(*(reinterpret_cast<const uint16_t*>(in)));
                break;
              case 1:
                *out = this->reverse(*in);
                break;
              default:
                BL_ERRORF("ERROR:  length of array should be > 0 and <= 8.  For longer, break it up or use the vector version of bitgroup_ops.");
                break;

            }

            return 0;
          }

          // bitgroup size is larger than the word type in bits.
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE>
          BITS_INLINE typename std::enable_if<(BITS > (sizeof(WORD_TYPE) << 3)), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) const {
            return u;
          }
          // bitgroup size is same size as the word type in bits.
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS == (sizeof(WORD_TYPE) << 3)), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) const {
            return u;
          }
          // this specialization will work for all WORD_TYPE sizes.
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE>
          BITS_INLINE typename std::enable_if<(BITS <= 8), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) const {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            WORD_TYPE v = u;

            // essentially a Duff's Device here, cascading to smaller swap sizes.
            // swap the bit_groups in logarithmic steps.  only execute steps that have swap size smaller than the WORD_TYPE, and only if bit group is <= the current swap size
            switch (sizeof(WORD_TYPE)) {
              case 8:
                v = _bliss_bswap_64(v);
                break;    // swap all 8 bytes
              case 4:
                v = _bliss_bswap_32(v);
                break;  // swap all 4 bytes
              case 2:
                v = _bliss_bswap_16(v);
                break;  // swap both 2 byte
              default:
                break;
            }

            // now swap if requested group size is 4, 2, or 1.  only execute if bit group is <= to current swap size.
            switch (BITS) {
              case 1:
                v = ((v >>  1) &  mask1) | ((v &  mask1) <<  1);  // swap 1 bits
              case 2:
                v = ((v >>  2) &  mask2) | ((v &  mask2) <<  2);  // swap 2 bits
              case 4:
                v = ((v >>  4) &  mask4) | ((v &  mask4) <<  4);  // swap nibbles
              default:
                break;
            }

            return v;
          }

          /// this specialization is for 8 byte words.   16 bit or 32 bit groups are allowed.
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS < (sizeof(WORD_TYPE) << 3)) && (sizeof(WORD_TYPE) == 8), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) const {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            WORD_TYPE v = u;

            // essentially a Duff's Device here, cascading to smaller swap sizes.
            // swap the bit_groups in logarithmic steps.  only execute steps that have swap size smaller than the WORD_TYPE, and only if bit group is <= the current swap size
            // BIT_GROUP_SIZE > 8:  16 or 32, depending on word size..  limited by word_type size
            switch (BITS) {
              // should be shifting by 32 but causes compiler warning if v is same size of BIT_GROUP_SIZE
              case 16:
                v = ((v >> 16) & mask16) | ((v & mask16) << 16);
              case 32:
                v = ((v >> 32) & mask32) | ((v & mask32) << 32);
              default:
                break;
            }

            return v;
          }

          /// this specialization is for 4 byte words.  only 16 bit group is allowed/needed here.
          /// separate by word_type size instead of bits so that we don't shift by 32 when word is 4 bytes.
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS < (sizeof(WORD_TYPE) << 3)) && (sizeof(WORD_TYPE) == 4), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) const {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            return ((u >> 16) & mask16) | ((u & mask16) << 16);
          }


          // for performance testing of the reverse_transform framework
          template <unsigned int BITS = BIT_GROUP_SIZE, typename WORD_TYPE>
          BITS_INLINE typename std::enable_if<(BITS < 8) && (BITS && (BITS - 1) == 0), WORD_TYPE>::type
          reverse_bits_in_byte(WORD_TYPE const &u) const {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            WORD_TYPE v = u;

            // now swap if requested group size is 4, 2, or 1.  only execute if bit group is <= to current swap size.
            switch (BITS) {
              case 1:
                v = ((v >>  1) &  mask1) | ((v &  mask1) <<  1);  // swap 1 bits
              case 2:
                v = ((v >>  2) &  mask2) | ((v &  mask2) <<  2);  // swap 2 bits
              case 4:
                v = ((v >>  4) &  mask4) | ((v &  mask4) <<  4);  // swap nibbles
              default:
                break;
            }

            return v;
          }

      };


      /// full template specialization for SWAR based bit reverse for 3 bits.  uses BSWAP when appropriate, else use bit shift.  this performs a SWAR style bit shuffling after a full bitwise reverse.)
      template <bool POW2>
      struct bitgroup_ops<3, BIT_REV_SWAR, POW2> {

          // masks.
          static constexpr uint64_t mask3lo = 0x9249249249249249;
          static constexpr uint64_t mask3mid = 0x2492492492492492;
          static constexpr uint64_t mask3hi = 0x4924924924924924;

          static constexpr uint64_t  mask_all = 0xFFFFFFFFFFFFFFFF;

          static constexpr unsigned int bitsPerGroup = 3;
          static constexpr unsigned char simd_type = BIT_REV_SWAR;

          bitgroup_ops<1, BIT_REV_SWAR, true> bit_rev_1;

          /// reverse function to reverse bits for data types are are not power of 2 number of bytes
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if (len == 0) return bit_offset; //throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            //if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be less than 8.  for larger, shift input instead.");
            assert(bit_offset < 8);

            // SHIFTS are NOT NECESSARY.
            // memcpy input is not needed.
            // NEED To zero the extra bytes.

            switch (len) {
              case 8:
                {
                  uint64_t v8 = *(reinterpret_cast<const uint64_t*>(in));
                  *(reinterpret_cast<uint64_t*>(out)) |= this->reverse(v8, bit_offset); // >> bit_offset;
                }
                break;
              case 7:
              case 6:
              case 5:
                {
                  uint64_t v7i = *(reinterpret_cast<const uint64_t*>(in));
                  v7i &= (mask_all >> (64 - (len << 3)));
                  uint64_t v7 = this->reverse(v7i, bit_offset);
                  out[0] |= *(reinterpret_cast<uint8_t*>(&v7) + 8 - len);
                  memcpy(out + 1, reinterpret_cast<uint8_t*>(&v7) + (9 - len), len - 2);
                  out[len - 1] |= *(reinterpret_cast<uint8_t*>(&v7) + 7);

                }
                break;
              case 4:
                {
                  uint32_t v4 = *(reinterpret_cast<const uint32_t*>(in));
                  *(reinterpret_cast<uint32_t*>(out)) |= this->reverse(v4, bit_offset); // >> bit_offset;
                }
                break;
              case 3:
                {
                  uint32_t v3i = *(reinterpret_cast<const uint32_t*>(in));
                  v3i &= (mask_all >> 40);
                  uint32_t v3 = this->reverse(v3i, bit_offset);
                  out[0] |= *(reinterpret_cast<uint8_t*>(&v3) + 1);
                  out[1] =  *(reinterpret_cast<uint8_t*>(&v3) + 2);
                  out[2] |= *(reinterpret_cast<uint8_t*>(&v3) + 3);

                }
                break;
              case 2:
                {
                  uint16_t v2 = *(reinterpret_cast<const uint16_t*>(in));
                  *(reinterpret_cast<uint16_t*>(out)) |= this->reverse(v2, bit_offset);
                }
                break;
              case 1:
                {
                  uint8_t v1 = *in;
                  *out |= this->reverse(v1, bit_offset);
                }
                break;
              default:
                BL_ERRORF("ERROR:  length of array should be > 0 and <= 8.  For longer, break it up or use the vector version of bitgroup_ops.");
                break;
            }

            return ((len << 3) - bit_offset) % 3;
          }

          // dispatcher
          template <typename WORD_TYPE>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & u, uint16_t const & bit_offset) const {
            switch (bit_offset % 3) {
              case 0: return reverse<0, WORD_TYPE>(::std::forward<WORD_TYPE const>(u)); break;
              case 1: return reverse<1, WORD_TYPE>(::std::forward<WORD_TYPE const>(u)); break;
              case 2: return reverse<2, WORD_TYPE>(::std::forward<WORD_TYPE const>(u)); break;
              default: return 0; break;
            }
          }

          /// reverse in groups of 3 bits.  NOTE: within the offset or remainder, the middle bit remains,
          /// while the high and low bits are zeroed during the reversal.  this makes OR'ing with adjacent
          /// entries simple.
          template <uint16_t offset = 0, typename WORD_TYPE>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const &u) const {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            // and swap the 1st and 3rd elements in the group.
            // we need 3 patterns for shifting based on the offset.
            // offset is at the LSB position of the ORIGINAL word

            WORD_TYPE v = u;
            // 2 shifts, 2 ors, 3 ands
            //switch (((sizeof(WORD_TYPE) << 3) - offset) % 3) {
            switch (offset % 3) {
              case 2:
                // rem == 2:                    mid, hi, lo
                v = ((v & mask3mid) >>  2) | (v & mask3lo ) | ((v & mask3hi ) << 2);
                break;
              case 1:
                // rem == 1:                    hi, lo, mid
                v = ((v & mask3lo ) >>  2) | (v & mask3hi ) | ((v & mask3mid) << 2);
                break;
              default:
                // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
                v = ((v & mask3hi ) >>  2) | (v & mask3mid) | ((v & mask3lo ) << 2);
                break;
            }

            //========================== finally reverse bits in groups of 1.
            return bit_rev_1.reverse(v);
            //============================ done reverse bits in groups of 1

          }

      };

#if defined(__SSSE3__)

      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
      /// DO THIS BECAUSE BITWISE SHIFT does not cross the epi64 boundary.
      template <>
      BITS_INLINE __m128i srli(__m128i const & val, uint16_t shift) {
        if (shift == 0) return val;
        if (shift >=128) return _mm_setzero_si128();

        // get the 9th byte into the 8th position
        // shift by 1 byte.  avoids a mask load + shuffle, and is friendly with epi16, epi32, and epi64 versions of shifts.
        __m128i tmp = _mm_srli_si128(val, 8);  //(bytes level shift only)


        // then left shift TMP to get the bits that crosses the 64 bit boundary.
        // shift the input value via epi64 by the number of shifts
        // then or together.
        if (shift < 64)
        	return _mm_or_si128(_mm_slli_epi64(tmp, 64 - shift), _mm_srli_epi64(val, shift));

        // 64 bit exactly
        if (shift == 64) return tmp;

//        if (shift > 64)
        	return _mm_srli_epi64(tmp, shift - 64);

      }

      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
      /// DO THIS BECAUSE BITWISE SHIFT does not cross the epi64 boundary.
      template <uint16_t SHIFT>
      BITS_INLINE __m128i srli(__m128i const & val) {
        if (SHIFT == 0) return val;
        if (SHIFT >=128) return _mm_setzero_si128();

        // get the 9th byte into the 8th position
        // shift by 1 byte.  avoids a mask load + shuffle, and is friendly with epi16, epi32, and epi64 versions of shifts.
        __m128i tmp = _mm_srli_si128(val, 8);  //(bytes level shift only)

        // then left shift TMP to get the bits that crosses the 64 bit boundary.
        // shift the input value via epi64 by the number of shifts
        // then or together.
        if (SHIFT < 64)
          return _mm_or_si128(_mm_slli_epi64(tmp, 64 - SHIFT), _mm_srli_epi64(val, SHIFT));
        // 64 bit exactly
        if (SHIFT == 64) return tmp;

//        if (SHIFT > 64)
          return _mm_srli_epi64(tmp, SHIFT - 64);

      }

      /// shift left by number of bits less than 8.
      template <>
      BITS_INLINE __m128i slli(__m128i const & val, uint16_t shift) {
        if (shift == 0) return val;
        if (shift >=128) return _mm_setzero_si128();

        // get the 9th byte into the 8th position
        // shift by 1 byte.  avoids a mask load during shuffle, and is friendly with epi16, epi32, and epi64 versions of shifts.
        __m128i tmp = _mm_slli_si128(val, 8);

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        if (shift < 64)
          return _mm_or_si128(_mm_srli_epi64(tmp, 64 - shift), _mm_slli_epi64(val, shift));

        if (shift > 64)
          return _mm_slli_epi64(tmp, shift - 64);

        // 64 bit exactly
        return tmp;
      }
      /// shift left by number of bits less than 8.
      template <uint16_t SHIFT>
      BITS_INLINE __m128i slli(__m128i const & val) {
      	if (SHIFT == 0) return val;
      	if (SHIFT >=128) return _mm_setzero_si128();

      	// get the 9th byte into the 8th position
        // shift by 1 byte.  avoids a mask load during shuffle, and is friendly with epi16, epi32, and epi64 versions of shifts.
        __m128i tmp = _mm_slli_si128(val, 8);

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        if (SHIFT < 64)
        	return _mm_or_si128(_mm_srli_epi64(tmp, 64 - SHIFT), _mm_slli_epi64(val, SHIFT));
        // 64 bit exactly
        if (SHIFT == 64) return tmp;

        // if (SHIFT > 64)
        return _mm_slli_epi64(tmp, SHIFT - 64);

      }


      template <>
      BITS_INLINE __m128i bit_not(__m128i const & u) {
        return _mm_xor_si128(u, _mm_cmpeq_epi8(u, u));  // no native negation operator, so use xor
      }

      template <typename DEST_WORD_TYPE, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<::std::is_same<DEST_WORD_TYPE, __m128i>::value, __m128i>::type
      loadu(WORD_TYPE const * u) {
    	  return _mm_loadu_si128(reinterpret_cast<__m128i const *>(u));
      }
      template <typename SRC_WORD_TYPE, typename WORD_TYPE,
       typename = typename ::std::enable_if<::std::is_same<SRC_WORD_TYPE, __m128i>::value>::type >
      BITS_INLINE void storeu(WORD_TYPE * u, __m128i const & val) {
    	  _mm_storeu_si128(reinterpret_cast<__m128i *>(u), val);
      }

      template <>
      BITS_INLINE __m128i bit_or(__m128i const & u, __m128i const & v) {
        return _mm_or_si128(u, v);
      }
      template <>
      BITS_INLINE __m128i bit_and(__m128i const & u, __m128i const & v) {
        return _mm_and_si128(u, v);
      }
      template <>
      BITS_INLINE __m128i bit_xor(__m128i const & u, __m128i const & v) {
        return _mm_xor_si128(u, v);
      }


      template <>
      BITS_INLINE __m128i zero() {
        __m128i tmp;

#if defined(__INTEL_COMPILER)
  #pragma warning push
  #pragma warning disable 592
#else  // last one is gcc, since everyone defines __GNUC__.  clang can use the same as well.
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuninitialized"
#endif
        tmp = _mm_xor_si128(tmp, tmp);
#if defined(__INTEL_COMPILER)
  #pragma warning pop
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic pop
#endif

        return tmp;
//        return _mm_setzero_si128();
      }
      template <>
      BITS_INLINE __m128i bit_max() {
        __m128i tmp;

#if defined(__INTEL_COMPILER)
  #pragma warning push
  #pragma warning disable 592
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuninitialized"
#endif
        tmp = _mm_cmpeq_epi8(tmp, tmp);
#if defined(__INTEL_COMPILER)
  #pragma warning pop
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic pop
#endif

        return tmp;
      }



      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 128bit)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < 128, "ERROR: BIT_GROUP_SIZE is greater than number of bits in __m128i");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");

          /// index for reversing the bytes in groups of varying size (8: 1 byte; 16: 2 bytes, etc.)
          const __m128i rev_idx8;
          const __m128i rev_idx16;
          /// mask for lower 4 bits
          const __m128i _mask_lo;
          /// lookup table for reversing bits in a byte in bit groups of 1
          const __m128i lut1_lo;
          const __m128i lut1_hi;
          /// lookup table for reversing bits in a byte in bit groups of 2
          const __m128i lut2_lo;
          const __m128i lut2_hi;

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_SSSE3;

          static ::std::string toString(__m128i const & v) {
            uint8_t BLISS_ALIGNED_ARRAY(tmp, 16, 16);
            _mm_store_si128((__m128i*)tmp, v);
            ::std::stringstream ss;
            for (int i = 15; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }

          bitgroup_ops() :
            rev_idx8(_mm_setr_epi8(0x0F,0x0E,0x0D,0x0C,0x0B,0x0A,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,0x00)),
            rev_idx16(_mm_setr_epi8(0x0E,0x0F,0x0C,0x0D,0x0A,0x0B,0x08,0x09,0x06,0x07,0x04,0x05,0x02,0x03,0x00,0x01)),
            _mask_lo(_mm_setr_epi8(0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F)),
            lut1_lo(_mm_setr_epi8(0x00,0x08,0x04,0x0c,0x02,0x0a,0x06,0x0e,0x01,0x09,0x05,0x0d,0x03,0x0b,0x07,0x0f)),
            lut1_hi(_mm_setr_epi8(0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0)),
            lut2_lo(_mm_setr_epi8(0x00,0x04,0x08,0x0c,0x01,0x05,0x09,0x0d,0x02,0x06,0x0a,0x0e,0x03,0x07,0x0b,0x0f)),
            lut2_hi(_mm_setr_epi8(0x00,0x40,0x80,0xc0,0x10,0x50,0x90,0xd0,0x20,0x60,0xa0,0xe0,0x30,0x70,0xb0,0xf0))
            {
          }

          /// reverse function to reverse bits for data types are are not __mm128i  (has to be bigger than at least half)
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if ((len << 3) < BIT_GROUP_SIZE) return bit_offset; // throw ::std::invalid_argument("ERROR reversing byte array: length is 0");
            //if (len > 16) throw std::invalid_argument("ERROR:  length of array should be <= 16.  For longer, use a different bitgroup_ops (AVX2) or break it up.");
            assert ( len <= 16);
//            if (((len << 3) % BIT_GROUP_SIZE) > 0)
//              throw ::std::invalid_argument("ERROR reversing byte array:  len in bits needs to be a multiple of BIT_GROUP_SIZE");
            assert(((len << 3) % BIT_GROUP_SIZE) == 0);

            //if (bit_offset > 0) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be 0 for SSSE3 bit reverse and power of 2 BIT_GROUP_SIZE.");
            assert(bit_offset == 0);

            // convert to __m128i, then call the __m128i version.  note that there aren't types that have 128 bits other than simd types

            if (len == 16) { // full 16 byte array and BIT_GROUP_SIZE is power of 2, so directly load and store.
            	_mm_storeu_si128((__m128i*)out, this->reverse( _mm_loadu_si128((__m128i*)in) ));
            } else { // not the full 16 bytes.  so need to copy out at the end to out array.
              // we can directly cast from __m128i to uint8_t

              // DO NOT NEED TO WORRY about reverse involving bits from outside the specified length, so directly loading from "in" is okay and no memset is needed.

              __m128i w = this->reverse(_mm_loadu_si128((__m128i*)in));

              // cast to uint8_t array for copy out
              uint8_t* tmp = reinterpret_cast<uint8_t*>(&w);

              // maskmoveu is very slow.  why?
              //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2
              memcpy(out, tmp + 16 - len, len);
            }
            return 0;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS >= 128), __m128i>::type
          reverse(__m128i const & u) const {
            static_assert(BIT_GROUP_SIZE == 128, "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }

          // groupsize =1 version:  4 loads, 3 shuffles, 3 logical ops, 1 shift op
          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS <= 8), __m128i>::type
          reverse(__m128i const & u) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            // if bit group is <= 8, then need to shuffle bytes, and then shuffle bits.
            __m128i v = _mm_shuffle_epi8(u, rev_idx8);

            //== now see if we need to shuffle bits.
            if (BIT_GROUP_SIZE == 8) return v;

            //== next get the lut indices.

            // lower 4 bits
            __m128i lo = _mm_and_si128(_mask_lo, v);                                                      // SSE2
            // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
            // shift not using si128 - that instruction shifts bytes not bits.
            __m128i hi = _mm_srli_epi16(_mm_andnot_si128(_mask_lo, v), 4);                                // SSE2

            switch (BIT_GROUP_SIZE) {
              case 1:
                //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                lo = _mm_shuffle_epi8(lut1_hi, lo);                        // SSSE3
                hi = _mm_shuffle_epi8(lut1_lo, hi);                        // SSSE3
                break;
              case 2:
                //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                lo = _mm_shuffle_epi8(lut2_hi, lo);                        // SSSE3
                hi = _mm_shuffle_epi8(lut2_lo, hi);                        // SSSE3
                break;
              case 4:
                lo = _mm_slli_epi16(lo, 4);  // shift lower 4 to upper
                break;
              default:
                break;

            }

            // recombine
            return _mm_or_si128(lo, hi);                                                                  // SSE2
          }

          // groupsize =1 version:  4 loads, 3 shuffles, 3 logical ops, 1 shift op
          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS < 128), __m128i>::type
          reverse(__m128i const & u) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            __m128i v;

            //== first reverse the bit groups via 1 shuffle operation
            switch (BIT_GROUP_SIZE) {
              case 16:   // if bit group is 16, then shuffle using epi8
                v = _mm_shuffle_epi8(u, rev_idx16);//TCP1
                break;                                // SSSE3
              case 32:
                v = _mm_shuffle_epi32(u, 0x1B);
                break; // original 11 10 01 00 (high to low) => 00 01 10 11 (high to low) == 0x1B                // SSE2
              case 64:
                v = _mm_shuffle_epi32(u, 0x4E);
                break;  // original 11 10 01 00 (high to low) => 01 00 11 10 (high to low) == 0x4E           // SSE2
              default:
                break;
            }

            return v;
          }

          // groupsize =1 version:  4 loads, 3 shuffles, 3 logical ops, 1 shift op
          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS < 8) && ((BITS & (BITS - 1)) == 0), __m128i>::type
          reverse_bits_in_byte(__m128i const & u) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

        	__m128i v = u;
            //== next get the lut indices.

            // lower 4 bits
            __m128i lo = _mm_and_si128(_mask_lo, v);                                                      // SSE2
            // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
            // shift not using si128 - that instruction shifts bytes not bits.
            __m128i hi = _mm_srli_epi16(_mm_andnot_si128(_mask_lo, v), 4);                                // SSE2

            switch (BIT_GROUP_SIZE) {
              case 1:
                //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                lo = _mm_shuffle_epi8(lut1_hi, lo);                        // SSSE3
                hi = _mm_shuffle_epi8(lut1_lo, hi);                        // SSSE3
                break;
              case 2:
                //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                lo = _mm_shuffle_epi8(lut2_hi, lo);                        // SSSE3
                hi = _mm_shuffle_epi8(lut2_lo, hi);                        // SSSE3
                break;
              case 4:
                lo = _mm_slli_epi16(lo, 4);  // shift lower 4 to upper
                break;
              default:
                break;

            }

            // recombine
            return _mm_or_si128(lo, hi);                                                                  // SSE2
          }
      };



      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 128bit)
      template <bool POW2>
      struct bitgroup_ops<3, BIT_REV_SSSE3, POW2> {

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = 3;
          static constexpr unsigned char simd_type = BIT_REV_SSSE3;

          const __m128i mask3lo;
          const __m128i mask3mid;
          const __m128i mask3hi;


          static constexpr uint8_t BLISS_ALIGNED_ARRAY(mask_all, 32, 16) = {0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
                                                                0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};

          bitgroup_ops<1, BIT_REV_SSSE3, true> bit_rev_1;

          bitgroup_ops() :
            mask3lo( _mm_setr_epi32(0x49249249, 0x92492492, 0x24924924, 0x49249249)),
            mask3mid(_mm_setr_epi32(0x92492492, 0x24924924, 0x49249249, 0x92492492)),
            mask3hi( _mm_setr_epi32(0x24924924, 0x49249249, 0x92492492, 0x24924924))
          {}

          static ::std::string toString(__m128i const & v) {
            uint8_t BLISS_ALIGNED_ARRAY(tmp, 16, 16);
            _mm_store_si128((__m128i*)tmp, v);
            ::std::stringstream ss;
            for (int i = 15; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }


          /// reverse function to reverse bits for data types are are not __mm128i  (has to be bigger than at least half)
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if (len == 0) return bit_offset; //throw ::std::invalid_argument("ERROR reversing byte array: length is 0");
            //if (len > 16) throw std::invalid_argument("ERROR:  length of array should be <= 16.  For longer, use a different bitgroup_ops (AVX2) or break it up.");
            assert (len <= 16);

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            //if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array SSSE3:  bit offset should be less than 8.  for larger, shift input instead.");
            assert (bit_offset < 8);

            // convert to __m128i, then call the __m128i version.  note that there aren't types that have 128 bits other than simd types

            // NOTE that reversed byte array will need to be OR'd together,

            __m128i v = _mm_loadu_si128((__m128i*)in);

            if (len < 16)
              // NOTE since reverse may involve bits outside of len, we have to zero the part outside of len.  only do if len is less than 16.
              v = _mm_and_si128(v, _mm_loadu_si128((__m128i*)(mask_all + 16 - len)));

            // ARE THE SHIFTS NECESSARY?  they end up zeroing out the offset, which could be accomplished via a mask.
            // note that the remainder portion is not zeroed, specifically, if rem == 2, then the center bit is not zeroed, but the side bits are zeroed.
            // center bit OR with itself does not change.  side bits from next iter ORed with zeroed side bits (this iter) returns the right results.
            // SO ARE THE ZEROING ALSO NECESSARY?

            // SHIFTS are NOT NECESSARY.

            // do the conversion
            v = this->reverse( v , bit_offset);

            // save the first and last byte, then copy the data back, finally OR the first and last byte back.
            uint8_t first = out[0], last = out[len-1];

            // copy the data
            if (len == 16) {
              _mm_storeu_si128((__m128i*)out, v);
            } else {  // has to do memcpy.
              // copy back to out. need to OR with existing output, so we always need to copy regardless of len.
              uint8_t * tmp = reinterpret_cast<uint8_t*>(&v);   // can reinterpret_cast directly.
              memcpy(out, tmp + (16 - len), len);
            }

            // or back the old values.
            out[0] |= first;
            out[len-1] |= last;

            // maskmoveu is very slow.  why?
			//_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2

            // return remainder.
            return ((len << 3) - bit_offset) % 3;

          }


          BITS_INLINE __m128i reverse(__m128i const & u, uint16_t bit_offset) const {
            switch (bit_offset % 3) {
              case 0: return reverse<0>(::std::forward<__m128i const>(u)); break;
              case 1: return reverse<1>(::std::forward<__m128i const>(u)); break;
              case 2: return reverse<2>(::std::forward<__m128i const>(u)); break;
              default: return _mm_setzero_si128();
                break;
            }
          }

          /// reverse in groups of 3 bits.  NOTE: within the offset or remainder, the middle bit remains,
          /// while the high and low bits are zeroed during the reversal.  this makes OR'ing with adjacent
          /// entries simple.
            template <uint16_t offset = 0>
            BITS_INLINE __m128i reverse(__m128i const & u ) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

              // get the individual bits.  3 loads, 3 ANDs.
              __m128i lo = _mm_and_si128(u, mask3lo);
              __m128i mid = _mm_and_si128(u, mask3mid);
              __m128i hi = _mm_and_si128(u, mask3hi);

              __m128i v;
              // next shift based on the rem bits.  2 cross boundary shifts = 2 byteshift, 4 shifts, 2 ORs. + 2 ORs.
              switch (offset % 3) {
                case 2:
                  // rem == 2:                    mid, hi, lo
                  v = _mm_or_si128(srli(mid, 2), _mm_or_si128(lo, slli( hi, 2) ) );
                  break;
                case 1:
                  // rem == 1:                    hi, lo, mid
                  v = _mm_or_si128(srli(lo, 2), _mm_or_si128(hi, slli( mid, 2) ) );
                  break;
                default:
                  // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
                  v = _mm_or_si128(srli(hi, 2), _mm_or_si128(mid, slli( lo, 2) ) );
                  break;
              }

              // alternatively,  for cross-boundary shift, do byte shift, and reverse shift (8-bit_offset), then or. with shifted (offset)->  each CBS is finished in 4 ops w/o load.


            //========================== first reverse bits in groups of 1.
            return bit_rev_1.reverse(v);
            //=========================== done reverse bits in groups of 1
          }


      };

      // need the DUMMY template parameter to correctly instantiate here.
      template <bool POW2>
      constexpr uint8_t bitgroup_ops<3, BIT_REV_SSSE3, POW2>::BLISS_ALIGNED_ARRAY(mask_all, 32, 16);



#endif

#if defined(__AVX2__)

      // ALSO try using memory...
      // also try using 64 bit shift

//      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
//      /// only support shift <= 1 byte at the moment, since alignr has to take an 8 bit immediate (i.e. const)
//      template <>
//      BITS_INLINE __m256i srli(__m256i const & val, uint16_t shift) {
//    	  if (shift == 0) return val;
//    	  if (shift >= 256) return _mm256_setzero_si256();
//
//    	// goal:  get next higher 64 bits to current position.
//    	  // conditionals:  shift byt 0 to 63, 64 to 128, 128 to 192, and 192 to 256
//    	  // code below can be modified to shift up to 64 bits, and with conditional, up to 128 bits?
//    	  // 64, 192: p+alignr is enough
//    	  // 128: permute is enough
//    	  // 0-63:  same as below - 5 instructions
//    	  // 65-127, 129-191: permute + 2 alignr, then or(slli, srli)  total 6 instructions
//    	  // 193-255: permute + alignr, then right shift. 3 instructions
//
//    	  //alternatives:  permute4x64 (all 256bits, but can't zero), s(l/r)li_si256 (128 bit lanes)
//    	 // alternative, do the 64 bit shift, then do remaining.
//
//        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
//        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
//        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
//        //                        to shift right, b = val, N = 1, lower lane of a should be upper lane of val, upper lane of a should be 0.  possibly 2 ops. permute latency is 3
//        __m256i hi = _mm256_permute2x128_si256(val, val, 0x83);  // 1000.0011 higher lane is 0, lower lane is higher lane of val
//        if (shift == 128) return hi;
//
//        uint16_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.
//        __m256i v1 = _mm256_alignr_epi8(hi, val, byte_shift); // alignr:  src1 low OR src2 low, then right shift bytes, put into output low.  same with high
//        uint16_t bit64_shift = shift & 0x3F;        // shift within the 64 bit block
//
//        if (bit64_shift == 0) return v1;  // exact byte alignment.  return it.
//        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
//        //    0   hi
//        //    hi  lo    => dest_lo =  hi lo shift
//        //                 dest_hi =   0 hi shift
//        if (byte_shift == 24) return _mm256_srli_epi64(v1, bit64_shift);
//
//        // else we need one more alignr
//        __m256i v2 = _mm256_alignr_epi8(hi, val, byte_shift + 8);
//
//        // then left shift to get the bits in the right place
//        // shift the input value via epi64 by the number of shifts
//        // then or together.
//        return _mm256_or_si256(_mm256_slli_epi64(v2, 64 - bit64_shift), _mm256_srli_epi64(v1, bit64_shift));
//      }

      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
      template <uint16_t SHIFT>
      BITS_INLINE __m256i srli(__m256i const & val) {
    	  if (SHIFT == 0) return val;
    	  if (SHIFT >= 256) return _mm256_setzero_si256();


        // val is  A B C D  (MSB to LSB)

        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
        //                        to shift left, a = val, N = 15, lower lane of b should be 0, higher lane of b should be lower lane of val.  possibly 2 ops, permute latency is 3.
        __m256i hi = _mm256_permute2x128_si256(val, val, 0x83);  // 1000.3 higher lane is 0, lower lane is higher lane of val latency = 3
        if (SHIFT == 128) return hi;   // hi is 0 0 A B

        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
        //     0  hi
        //    hi  lo    => dest_lo =  hi lo shift
        //                 dest_hi =   0 hi shift
        //uint16_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.

        __m256i v1 = (SHIFT < 128) ? _mm256_alignr_epi8(hi, val, 8) :  // shift by 8
        //  else shift > 128
//                   _mm256_permute4x64_epi64(hi, 0xFD); // permute and shift hi right more.  3 3 3 1 -> 0xFD   latency = 3
                   //_mm256_alignr_epi8(zero, hi, 8);    // latency of 1  this is correct, but for consistency, use below
            _mm256_srli_si256(hi, 8);  // shift high to right by 8 bytes more, total 24 bytes.  latency 1.

    // if shift < 128, v1 is 0 A B C.  else v1 is 0 0 0 A

        constexpr uint8_t bit64_shift = SHIFT & 0x3F;        // shift within the 64 bit block

        if (bit64_shift == 0) return v1;  // already got cases at 0, 128, and 256.  this catches 64 and 192

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        __m256i lv, rv;
        switch (SHIFT >> 6) {
      case 0:
        rv = v1; lv = val; break;  // shift is < 64.  shift ABCD left and BCD0 right
      case 1:
        rv = hi; lv = v1; break;  // shift is < 128.  shift BCD0 left and CD00 right
      case 2:
        rv = v1; lv = hi; break;  // shift is < 192.  shift CD00 left and D000 right
      case 3:
        return _mm256_srli_epi64(v1, bit64_shift);  // shift is > 192, shift D000 left
      default:
        break;
        };

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        return _mm256_or_si256(_mm256_slli_epi64(rv, 64 - bit64_shift), _mm256_srli_epi64(lv, bit64_shift));


      }



//      /// shift left by number of bits less than 8.
//      template <>
//      BITS_INLINE __m256i slli(__m256i const & val, uint16_t shift) {
//    	  if (shift == 0) return val;
//    	  if (shift >= 256) return _mm256_setzero_si256();
//
//    	  // val is  A B C D
//
//        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
//        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
//        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
//        //                        to shift left, a = val, N = 15, lower lane of b should be 0, higher lane of b should be lower lane of val.  possibly 2 ops, permute latency is 3.
//        __m256i lo = _mm256_permute2x128_si256(val, val, 0x08);  // lower lane is 0, higher lane is lower lane of val
//        if (shift == 128) return lo;   // lo is C D 0 0
//
//        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
//        //    hi  lo
//        //    lo   0    => dest_lo =  hi lo shift
//        //                 dest_hi =   0 hi shift
//        //uint8_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.
//
//        // TODO - left shift is trickier.
//		__m256i v1 = (shift < 128) ? _mm256_alignr_epi8(val, lo, 8) :  // permute essentially shifted by 16
//				//  else shift > 128
//									 _mm256_permute4x64_epi64(lo, 0x80); // permute and shift lo left more.  2 0 0 0 -> 0x80
//		// if shift < 128, v1 is B C D 0.  else v1 is D 0 0 0
//
//        uint8_t bit64_shift = shift & 0x3F;        // shift within the 64 bit block
//
//        if (bit64_shift == 0) return v1;
//
//        // then left shift to get the bits in the right place
//        // shift the input value via epi64 by the number of shifts
//        // then or together.
//        __m256i lv, rv;
//        switch (shift >> 6) {
//			case 0:
//				rv = val; lv = v1; break;  // shift is < 64.  shift ABCD left and BCD0 right
//			case 1:
//				rv = v1; lv = lo; break;  // shift is < 128.  shift BCD0 left and CD00 right
//			case 2:
//				rv = lo; lv = v1; break;  // shift is < 192.  shift CD00 left and D000 right
//			case 3:
//				return _mm256_slli_epi64(v1, bit64_shift);  // shift is > 192, shift D000 left
//			default:
//				break;
//        };
//
//        // then left shift to get the bits in the right place
//        // shift the input value via epi64 by the number of shifts
//        // then or together.
//        return _mm256_or_si256(_mm256_srli_epi64(lv, 64 - bit64_shift), _mm256_slli_epi64(rv, bit64_shift));
//      }


      /// shift left by number of bits less than 8.
      template <uint16_t SHIFT>
      BITS_INLINE __m256i slli(__m256i const & val) {
    	  if (SHIFT == 0) return val;
    	  if (SHIFT >= 256) return _mm256_setzero_si256();

    	  // val is  A B C D (MSB to LSB)

        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
        //                        to shift left, a = val, N = 15, lower lane of b should be 0, higher lane of b should be lower lane of val.  possibly 2 ops, permute latency is 3.
        __m256i lo = _mm256_permute2x128_si256(val, val, 0x08);  // 0.1000 lower lane is 0, higher lane is lower lane of val  latency 3
        if (SHIFT == 128) return lo;   // lo is C D 0 0

        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
        //    hi  lo
        //    lo   0    => dest_lo =  hi lo shift
        //                 dest_hi =   0 hi shift
        //uint8_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.

		__m256i v1 = (SHIFT < 128) ? _mm256_alignr_epi8(val, lo, 8) :  // left shift val 8 bytes
				//  else shift > 128
		    //_mm256_permute4x64_epi64(lo, 0x80); // left shift lo 8 bytes.  2 0 0 0 -> 0x80 (3rd element move to 4th pos.  latency 3
		    //_mm256_alignr_epi8(lo, zero, 8);   // left shift by 24 bytes.  combine lo and zero then shift by 8.  latency 1
		    _mm256_slli_si256(lo, 8);  // left shift by 8 bytes within the lane to get effect of shift by 24 bytes total.  latency 1
		// if shift < 128, v1 is B C D 0.  else v1 is D 0 0 0

        constexpr uint8_t bit64_shift = SHIFT & 0x3F;        // shift within the 64 bit block

        if (bit64_shift == 0) return v1;  // already got cases at 0, 128, and 256.  this catches 64 and 192

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        __m256i lv, rv;
        switch (SHIFT >> 6) {
			case 0:
				rv = val; lv = v1; break;  // shift is < 64.  shift ABCD left and BCD0 right
			case 1:
				rv = v1; lv = lo; break;  // shift is < 128.  shift BCD0 left and CD00 right
			case 2:
				rv = lo; lv = v1; break;  // shift is < 192.  shift CD00 left and D000 right
			case 3:
				return _mm256_slli_epi64(v1, bit64_shift);  // shift is > 192, shift D000 left
			default:
				break;
        };

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        return _mm256_or_si256(_mm256_srli_epi64(lv, 64 - bit64_shift), _mm256_slli_epi64(rv, bit64_shift));
      }

      template <>
      BITS_INLINE __m256i bit_not(__m256i const & u) {
        return _mm256_xor_si256(u, _mm256_cmpeq_epi8(u, u));  // no native negation operator
      }

      template <typename DEST_WORD_TYPE, typename WORD_TYPE>
      BITS_INLINE typename ::std::enable_if<::std::is_same<DEST_WORD_TYPE, __m256i>::value, __m256i>::type
      loadu(WORD_TYPE const * u) {
    	  return _mm256_loadu_si256(reinterpret_cast<__m256i const *>(u));
      }
      template <typename SRC_WORD_TYPE, typename WORD_TYPE,
       typename = typename ::std::enable_if<::std::is_same<SRC_WORD_TYPE, __m256i>::value>::type >
      BITS_INLINE void storeu(WORD_TYPE * u, __m256i const & val) {
    	  _mm256_storeu_si256(reinterpret_cast<__m256i *>(u), val);
      }

      template <>
      BITS_INLINE __m256i bit_or(__m256i const & u, __m256i const & v) {
        return _mm256_or_si256(u, v);
      }
      template <>
      BITS_INLINE __m256i bit_and(__m256i const & u, __m256i const & v) {
        return _mm256_and_si256(u, v);
      }
      template <>
      BITS_INLINE __m256i bit_xor(__m256i const & u, __m256i const & v) {
        return _mm256_xor_si256(u, v);
      }
      template <>
      BITS_INLINE __m256i zero() {
//        return _mm256_setzero_si256();
        __m256i tmp;
#if defined(__INTEL_COMPILER)
  #pragma warning push
  #pragma warning disable 592
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuninitialized"
#endif
        tmp = _mm256_xor_si256(tmp, tmp);
#if defined(__INTEL_COMPILER)
  #pragma warning pop
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic pop
#endif
        return tmp;
      }
      template <>
      BITS_INLINE __m256i bit_max() {
        __m256i tmp;

#if defined(__INTEL_COMPILER)
  #pragma warning push
  #pragma warning disable 592
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wuninitialized"
#endif
        tmp = _mm256_cmpeq_epi8(tmp, tmp);
#if defined(__INTEL_COMPILER)
  #pragma warning pop
#else  // last one is gcc, since everyone defines __GNUC__
  #pragma GCC diagnostic pop
#endif

        return tmp;

      }

      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 256bit)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < 256, "ERROR: BIT_GROUP_SIZE is greater than number of bits in __m256i");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");

          static const __m256i rev_idx_lane8;
          static const __m256i rev_idx_lane16;
          static const __m256i rev_idx32;
          static const __m256i _mask_lo;
          static const __m256i lut1_lo;
          static const __m256i lut1_hi;
          static const __m256i lut2_lo;
          static const __m256i lut2_hi;
          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_AVX2;

          bitgroup_ops() {}

          static ::std::string toString(__m256i const & v) {
            uint8_t BLISS_ALIGNED_ARRAY(tmp, 32, 32);
            _mm256_store_si256((__m256i*)tmp, v);																// AVX
            ::std::stringstream ss;
            for (int i = 31; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }


          /// reverse function to reverse bits for data types are are not __mm256i   (has to be bigger than at least half)
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if ((len << 3) < BIT_GROUP_SIZE) return bit_offset; // throw ::std::invalid_argument("ERROR reversing byte array: length is 0");
            //if (len > 32) throw std::invalid_argument("ERROR:  length of array should be <= 32.  For longer, break it up.");
            assert(len <= 32);

//            if (((len << 3) % BIT_GROUP_SIZE) > 0)
//              throw ::std::invalid_argument("ERROR reversing byte array:  len in bits needs to be a multiple of BIT_GROUP_SIZE");
            assert(((len << 3) % BIT_GROUP_SIZE) == 0);

            //if (bit_offset > 0) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be 0 for AVX2 bit reverse and power of 2 BIT_GROUP_SIZE.");
            assert(bit_offset == 0);

            // convert to __m256i, then call the __m256i version.  note that there aren't types that have 256 bits other than simd types

            if (len == 32) {  // full 32 byte array and BIT_GROUP_SIZE is power of 2, so directly load and store.
              // copy in, do reverse, copy back out
              // NOTE: aligned load does not seem to work correctly - SEGV.  use unaligned load.
              _mm256_storeu_si256((__m256i*)out, this->reverse( _mm256_loadu_si256((__m256i*)in) ));                    // AVX

            } else {  // not the full 32 bytes.  so need to copy out at the end to out array.
              // we can directly cast from __m256i to uint8_t

              // DO NOT NEED TO WORRY about reverse involving bits from outside the specified length, so directly loading from "in" is okay and no memset is needed.

              // copy in, do reverse, copy back out
              // NOTE: aligned load does not seem to work correctly - SEGV.  use unaligned load.
              __m256i w = this->reverse( _mm256_loadu_si256((__m256i*)in) );

              // cast to uint8_t array for copy out
              uint8_t *tmp = reinterpret_cast<uint8_t*>(&w);

              // maskmoveu is very slow.  why?
              //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2
              memcpy(out, tmp + 32 - len, len);
            }
            return 0;
          }



          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS >= 256), __m256i>::type
          reverse(__m256i const & u) const {
            static_assert(BIT_GROUP_SIZE == 256, "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS <= 8), __m256i>::type
          reverse(__m256i const & u) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            __m256i v = u;

            // if bit group is <= 8, then need to shuffle bytes, and then shuffle bits.
            v = _mm256_shuffle_epi8(v, rev_idx_lane8);     // reverse 8 in each of 128 bit lane          // AVX2
            v = _mm256_permute2x128_si256(v, v, 0x03);		// then swap the lanes                           // AVX2.  latency = 3

            //== now see if we need to shuffle bits.
            if (BITS == 8) return v;
			//== next get the lut indices.

			// lower 4 bits
			__m256i lo = _mm256_and_si256(_mask_lo, v); // no and_si256.  had to cast to ps first.           // AVX2
			// upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
			// shift not using si128 - that instruction shifts bytes not bits.
			__m256i hi = _mm256_srli_epi16(_mm256_andnot_si256(_mask_lo, v), 4);                                // AVX2

			switch (BIT_GROUP_SIZE) {
			case 1:
			  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
			  lo = _mm256_shuffle_epi8(lut1_hi, lo);                        // AVX2
			  hi = _mm256_shuffle_epi8(lut1_lo, hi);                        // AVX2
			  break;
			case 2:
			  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
			  lo = _mm256_shuffle_epi8(lut2_hi, lo);                        // AVX2
			  hi = _mm256_shuffle_epi8(lut2_lo, hi);                        // AVX2
			  break;
			case 4:
			  lo = _mm256_slli_epi16(lo, 4);  // shift lower 4 to upper										//AVX2
			  break;
			default:
			  break;
			}

			// recombine
			return _mm256_or_si256(lo, hi);                                                                  // AVX2

          }


          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS < 256), __m256i>::type
          reverse(__m256i const & u) const {


            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            __m256i v = u;

            //== first reverse the bit groups via 1 shuffle operation
            switch (BIT_GROUP_SIZE) {
              case 16:
                v = _mm256_shuffle_epi8(v, rev_idx_lane16);     // reverse 16 in each of 128 bit lane         // AVX2
                v = _mm256_permute2x128_si256(v, v, 0x03);							    // then swap the lanes							 // AVX2  latency = 3
                break;
              case 32:
                v = _mm256_permutevar8x32_epi32(v, rev_idx32);                         // AVX2 latency = 3
                break; // original 11 10 01 00 (high to low) => 00 01 10 11 (high to low) == 0x1B          // AVX2
              case 64:
                v = _mm256_permute4x64_epi64(v, 0x1B);                                                              // AVX2 latency = 3
                break; // original 11 10 01 00 (high to low) => 00 01 10 11 (high to low) == 0x1B         // AVX2
              case 128:
                v = _mm256_permute2x128_si256(v, v, 0x03);                                                         // AVX2 latency = 3
                break; //  low half of first -> hi, high half of second -> low  (coding is 0 1 2 3 => a_lo, a_hi, b_lo, b_hi)  // AVX2
              default:
                break;
            }

            return v;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS < 8) && (BITS && (BITS - 1) == 0), __m256i>::type
          reverse_bits_in_byte(__m256i const & u) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            __m256i v = u;
			//== next get the lut indices.

			// lower 4 bits
			__m256i lo = _mm256_and_si256(_mask_lo, v); // no and_si256.  had to cast to ps first.           // AVX2
			// upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
			// shift not using si128 - that instruction shifts bytes not bits.
			__m256i hi = _mm256_srli_epi16(_mm256_andnot_si256(_mask_lo, v), 4);                                // AVX2

			switch (BIT_GROUP_SIZE) {
			case 1:
			  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
			  lo = _mm256_shuffle_epi8(lut1_hi, lo);                        // AVX2
			  hi = _mm256_shuffle_epi8(lut1_lo, hi);                        // AVX2
			  break;
			case 2:
			  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
			  lo = _mm256_shuffle_epi8(lut2_hi, lo);                        // AVX2
			  hi = _mm256_shuffle_epi8(lut2_lo, hi);                        // AVX2
			  break;
			case 4:
			  lo = _mm256_slli_epi16(lo, 4);  // shift lower 4 to upper										//AVX2
			  break;
			default:
			  break;
			}

			// recombine
			return _mm256_or_si256(lo, hi);                                                                  // AVX2

          }


      };
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::rev_idx_lane8 =  _mm256_setr_epi32(0x0C0D0E0F, 0x08090A0B, 0x04050607, 0x00010203, 0x0C0D0E0F, 0x08090A0B, 0x04050607, 0x00010203);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::rev_idx_lane16 = _mm256_setr_epi32(0x0D0C0F0E, 0x09080B0A, 0x05040706, 0x01000302, 0x0D0C0F0E, 0x09080B0A, 0x05040706, 0x01000302);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::rev_idx32 =      _mm256_setr_epi32(0x00000007, 0x00000006, 0x00000005, 0x00000004, 0x00000003, 0x00000002, 0x00000001, 0x00000000);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::_mask_lo =       _mm256_setr_epi32(0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut1_lo = _mm256_setr_epi32(0x0c040800, 0x0e060a02, 0x0d050901, 0x0f070b03, 0x0c040800, 0x0e060a02, 0x0d050901, 0x0f070b03);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut1_hi = _mm256_setr_epi32(0xc0408000, 0xe060a020, 0xd0509010, 0xf070b030, 0xc0408000, 0xe060a020, 0xd0509010, 0xf070b030);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut2_lo = _mm256_setr_epi32(0x0c080400, 0x0d090501, 0x0e0a0602, 0x0f0b0703, 0x0c080400, 0x0d090501, 0x0e0a0602, 0x0f0b0703);
template <unsigned int BIT_GROUP_SIZE, bool POW2> const __m256i bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut2_hi = _mm256_setr_epi32(0xc0804000, 0xd0905010, 0xe0a06020, 0xf0b07030, 0xc0804000, 0xd0905010, 0xe0a06020, 0xf0b07030);

      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 128bit)
      template <bool POW2>
      struct bitgroup_ops<3, BIT_REV_AVX2, POW2> { 

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = 3;
          static constexpr unsigned char simd_type = BIT_REV_AVX2;

          static constexpr uint8_t BLISS_ALIGNED_ARRAY(mask_all, 64, 32) = {0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
                                                                0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
                                                                0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
                                                                0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};

          bitgroup_ops<1, BIT_REV_AVX2, true> bit_rev_1;

          static const __m256i mask3lo;
          static const __m256i mask3mid;
          static const __m256i mask3hi;


          bitgroup_ops() {}


          static ::std::string toString(__m256i const & v) {
            uint8_t BLISS_ALIGNED_ARRAY(tmp, 32, 32);
            _mm256_store_si256((__m256i*)tmp, v);
            ::std::stringstream ss;
            for (int i = 31; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }



          /// reverse function to reverse bits for data types are are not __mm128i  (has to be bigger than at least half)
          BITS_INLINE uint16_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint16_t const bit_offset = 0) const {
            if (len == 0) return bit_offset;//throw ::std::invalid_argument("ERROR reversing byte array: length is 0");
            //if (len > 32) throw std::invalid_argument("ERROR:  length of array should be <= 32.  For longer, break it up.");
            assert(len <= 32);

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            //if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array AVX2:  bit offset should be less than 8.  for larger, shift input instead.");
            assert(bit_offset < 8);

            // convert to __m256i, then call the __m256i version.  note that there aren't types that have 256 bits other than simd types

            // NOTE that reversed byte array will need to be OR'd together, 

            __m256i v = _mm256_loadu_si256((__m256i*)in);
            
            if (len < 32)  // only need to do this if len is < 32
              // NOTE since reverse may involve bits outside of len, we have to zero the part outside of len.
              v = _mm256_and_si256(v, _mm256_loadu_si256((__m256i*)(mask_all + 32 - len)));

            // do the conversion
            // ARE THE SHIFTS NECESSARY?  they end up zeroing out the offset, which could be accomplished via a mask.
            // note that the remainder portion is not zeroed, specifically, if rem == 2, then the center bit is not zeroed, but the side bits are zeroed.
            // center bit OR with itself does not change.  side bits from next iter ORed with zeroed side bits (this iter) returns the right results.
            // SO ARE THE ZEROING ALSO NECESSARY?

            // SHIFTS are NOT NECESSARY.
            v = this->reverse( v, bit_offset );


            // save the first and last byte, then copy the data back, finally OR the first and last byte back.
            uint8_t first = out[0], last = out[len-1];

            // copy the data
            if (len == 32) {
              _mm256_storeu_si256((__m256i*)out, v);
            } else {  // has to do memcpy.
              // copy back to out. need to OR with existing output, so we always need to copy regardless of len.
              uint8_t * tmp = reinterpret_cast<uint8_t*>(&v);   // can reinterpret_cast directly.
              memcpy(out, tmp + (32 - len), len);
            }

            // or back the old values.
            out[0] |= first;
            out[len-1] |= last;

            // maskmoveu is very slow.  why?
			     //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2

            // return remainder.
            return ((len << 3) - bit_offset) % 3;
          }

          BITS_INLINE __m256i reverse(__m256i const & u, uint16_t bit_offset) const {
            switch (bit_offset % 3) {
              case 0: return reverse<0>(::std::forward<__m256i const>(u)); break;
              case 1: return reverse<1>(::std::forward<__m256i const>(u)); break;
              case 2: return reverse<2>(::std::forward<__m256i const>(u)); break;
              default: return _mm256_setzero_si256();
                break;
            }
          }

          /// reverse in groups of 3 bits.  NOTE: within the offset or remainder, the middle bit remains,
          /// while the high and low bits are zeroed during the reversal.  this makes OR'ing with adjacent
          /// entries simple.
          template <uint16_t offset = 0>
          BITS_INLINE __m256i reverse(__m256i const & u) const {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            // get the individual bits.  3 loads, 3 ORs.
            __m256i lo = _mm256_and_si256(u, mask3lo);
            __m256i mid = _mm256_and_si256(u, mask3mid);
            __m256i hi = _mm256_and_si256(u, mask3hi);

            __m256i v;
            // next shift based on the rem bits.  2 cross boundary shifts = 2 byteshift, 4 shifts, 2 ORs. + 2 ORs.
            switch (offset % 3) {
            case 2:
              // rem == 2:                    mid, hi, lo
              v = _mm256_or_si256(srli<2>(mid), _mm256_or_si256(lo, slli<2>(hi)) );
              break;
            case 1:
              // rem == 1:                    hi, lo, mid
                v = _mm256_or_si256(srli<2>(lo), _mm256_or_si256(hi, slli<2>(mid)) );
              break;
            default:
              // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
                v = _mm256_or_si256(srli<2>(hi), _mm256_or_si256(mid, slli<2>(lo)) );
              break;
            }
            // alternatively,  for cross-boundary shift, do byte shift, and reverse shift (8-bit_offset), then or. with shifted (offset)->  each CBS is finished in 4 ops w/o load.

            //========================== first reverse bits in groups of 1.
            return bit_rev_1.reverse(v);
            //=========================== done reverse bits in groups of 1.

          }


      };
      template <bool POW2> const __m256i bitgroup_ops<3, BIT_REV_AVX2, POW2>::mask3lo  = 
        _mm256_setr_epi32(0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492);
      template <bool POW2> const __m256i bitgroup_ops<3, BIT_REV_AVX2, POW2>::mask3mid = 
        _mm256_setr_epi32(0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924);
      template <bool POW2> const __m256i bitgroup_ops<3, BIT_REV_AVX2, POW2>::mask3hi  =
        _mm256_setr_epi32(0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249);
      template <bool POW2> constexpr uint8_t bitgroup_ops<3, BIT_REV_AVX2, POW2>::BLISS_ALIGNED_ARRAY(mask_all, 64, 32);



#endif

      /**
       * @brief dispatcher for different bitgroup_ops specializations.  determine the best SIMD_TYPE given the WORD_TYPE, NWORDS, and BIT_GROUP_SIZE
       * @details  requirements:  word type should be unsigned integral,  BIT_GROUP_SIZE should not be more than nbits.
       *                  |       nbits      |   BIT_GROUP_SIZE  |   compiler flag  |
       *            SEQ   |       any        |        any        |      none        |
       *            SWAR  |       any        |     > 8; pow2, 3  |      none        |
       *            BSWAP |       any        |    <= 8; pow2, 3  |      none        |
       *            SSSE3 |     > 64         |        pow2, 3    |   SSSE3, AVX2    |
       *            AVX2  |     > 128        |        pow2, 3    |      AVX2        |
       *
       *    if the requirement is not satisfied, then it should fallback to next lowest.
       */

      // no arbitrary BIT_GROUP_SIZE, therefore no SEQ bitgroup_ops supported.

      // 2 sets of functions:  fixed sized array, and arbitrary length (data via pointer)

      template <typename WORD_TYPE, size_t len>
      void print(WORD_TYPE const (&words)[len]) {
        for (int i = len - 1; i >= 0; --i) {
          std::cout << std::hex << std::setfill('0') << std::setw(sizeof(WORD_TYPE) * 2) << words[i] << " ";
        }
      }
      template <typename WORD_TYPE>
      void print(WORD_TYPE const & word) {
        uint64_t const * ptr = reinterpret_cast<uint64_t const *>(&word);
        // uint64_t is 8 bytes, or shift by 3.
        for (int i = ((sizeof(WORD_TYPE) + sizeof(uint64_t) - 1) >> 3) - 1; i >= 0; --i) {
          std::cout << std::hex << std::setfill('0') << std::setw(std::min(sizeof(WORD_TYPE), sizeof(uint64_t)) * 2) << ptr[i] << " ";
        }
      }


      //============= generic transform operations ==================
      //  supports unary and binary transforms.
      //  4 possible variants:  a. no shift and no reverse
      //                        b. shift, no reverse
      //                        c. reverse, no shift
      //                        d. shift and reverse
      //  bitwise operation falls under a
      //  left and right shift operations falls under b
      //  a is a special case of b, so use b signature.
      //  c is a special case of d. so use d signature
      //  c and d are only for reversing bits/bit groups,
      //   note that the functor is responsible for reversing
      //    the bit groups within a machine word.  shift exists here because
      //    of padding bits.
      // TODO: we will ignore the case when we reverse and do true shift, for now.

      //========= unary transform - no shift======================

      // WORD orders bit bitwise transforms should be same between input and output
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)) , int>::type = 1>
      BITS_INLINE void bit_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
//    	  printf("0 bit transform exact. %lu x %lu byte words, machine word size %lu\n", len, sizeof(WORD_TYPE), sizeof(MachineWord));
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          storeu<MachineWord>(out, op(loadu<MachineWord>(in)));
      }
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE void bit_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

          //printf("0 bit transform smaller. %lu x %lu byte words, machine word size %lu\n", len, sizeof(WORD_TYPE), sizeof(MachineWord));

          MachineWord x = loadu<MachineWord>(in);
          x = op(x);
          store_part<(sizeof(WORD_TYPE) * len), 0>(out, x);

      }
      // no shift.
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE void bit_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
//    	  printf("0 bit transform by large\n");
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);

        MachineWord x = loadu<MachineWord>(u);  // faster than memcpy
        MachineWord s = op(x);

        // handle remainders first
        if (rem > 0) {
        	u += rem;

        	x = loadu<MachineWord>(u);

        	storeu<MachineWord>(w, s);
        	w += rem;

        	s = op(x);
        }
    	u += sizeof(MachineWord);

        // no overlap
        for (unsigned int i = 1; i < nMachWord; ++i) {
        	// load the next
        	x = loadu<MachineWord>(u);
            u += sizeof(MachineWord);

            // write the current
            storeu<MachineWord>(w, s);
            w += sizeof(MachineWord);

            // compute next
            s = op(x);
        }
        storeu<MachineWord>(w, s);
      }

      //=================== binary transform - no shift. ===========================

      // WORD orders bit bitwise transforms should be same between input and output, just shifted
      // positive shift increase value, so shift left.  negative shift shifts right. 0 shift returns original value, not sure if it's no-op.
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE void bit_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len], OP const & op) {

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          storeu<MachineWord>(out, op(loadu<MachineWord>(lhs), loadu<MachineWord>(rhs)));
      }
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)) , int>::type = 1>
      BITS_INLINE void bit_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len], OP const & op) {

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

          MachineWord x = op(loadu<MachineWord>(lhs), loadu<MachineWord>(rhs));
          store_part<bytes, 0>(out, x);
      }
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE void bit_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len], OP const & op) {

        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;


        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(lhs);
        uint8_t const * v = reinterpret_cast<uint8_t const *>(rhs);

        MachineWord x = loadu<MachineWord>(u);  // faster than memcpy
        MachineWord y = loadu<MachineWord>(v);  // faster than memcpy
        MachineWord s = op(x, y);

        // handle remainders first
        if (rem > 0) {
        	u += rem;
        	v += rem;

        	x = loadu<MachineWord>(u);
        	y = loadu<MachineWord>(v);

        	storeu<MachineWord>(w, s);
        	w += rem;

        	s = op(x, y);
        }
    	u += sizeof(MachineWord);
    	v += sizeof(MachineWord);

        // no overlap
        for (unsigned int i = 1; i < nMachWord; ++i) {
        	// load the next
        	x = loadu<MachineWord>(u);
        	y = loadu<MachineWord>(v);
            u += sizeof(MachineWord);
            v += sizeof(MachineWord);

            // write the current
            storeu<MachineWord>(w, s);
            w += sizeof(MachineWord);

            // compute next
            s = op(x, y);
        }
        storeu<MachineWord>(w, s);
      }

      // ========== unary transform with shift. =============================.
      // op being an identity does not add overhead.
      // WORD orders bit bitwise transforms should be same between input and output, just shifted
      // positive shift increase value, so shift left.  negative shift shifts right. 0 shift returns original value, not sure if it's no-op.

      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(BIT_SHIFT == 0), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
//    	  printf("0 shift transform by %d\n", BIT_SHIFT);
    	  //  		  std::cout << "in x: "; print(in); std::cout << std::endl;
    	  bit_transform<MAX_SIMD_TYPE, WORD_TYPE, len>(out, in, op);
//          std::cout << "shifted x: "; print(out); std::cout << std::endl;
      }

      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
          (BIT_SHIFT > 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert(BIT_SHIFT <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
//    	  printf("left shift transform, exact, by %d\n", BIT_SHIFT);
    	  using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          storeu<MachineWord>(out, ::bliss::utils::bit_ops::slli<BIT_SHIFT>(op(loadu<MachineWord>(in))));
      }
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
          (BIT_SHIFT < 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert((-BIT_SHIFT) <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
    	  using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
//    	  printf("right shift transform, exact, by %d, machineword bytes %lu\n", BIT_SHIFT, sizeof(MachineWord));
    	  MachineWord x = loadu<MachineWord>(in);
//    	  std::cout << " in: "; print(x); std::cout << std::endl;
    	  x = ::bliss::utils::bit_ops::srli<(0 - BIT_SHIFT)>(op(x));
          storeu<MachineWord>(out, x);
      }
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
                  (BIT_SHIFT > 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert(BIT_SHIFT <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
//    	  printf("left shift transform, small, by %d\n", BIT_SHIFT);
        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        // left shift.  okay.
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          // MachineWord is bigger than input, so there are extra bytes.  if reading from in, then right shift would result in junk on the left.
          // instead, for right shift, read so that the input are at high bytes of machine word.

          MachineWord x = ::bliss::utils::bit_ops::slli<BIT_SHIFT>(op(loadu<MachineWord>(in)));
          store_part<bytes, 0>(out, x);
      }
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
                  (BIT_SHIFT < 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert((-BIT_SHIFT) <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
//    	  printf("right shift transform, small, by %d\n", BIT_SHIFT);
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
    	  constexpr size_t bytes = sizeof(WORD_TYPE) * len;
          // right shift, load to upper bits.
          constexpr size_t offset = sizeof(MachineWord) - bytes;
          // MachineWord is bigger than input, so there are extra bytes.  if reading from in, then right shift would result in junk on the left.
          // instead, for right shift, read so that the input are at high bytes of machine word.

          MachineWord x =
        		  ::bliss::utils::bit_ops::srli<(-BIT_SHIFT)>(op(loadu<MachineWord>(reinterpret_cast<uint8_t const *>(in) - offset)));
          store_part<bytes, offset>(out, x);
      }

      // left shift byte aligned.  safe for in and out to be same.
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
          ((BIT_SHIFT & 0x7) == 0) && (BIT_SHIFT > 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert(BIT_SHIFT <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
//    	  printf("left shift transform, large, bytealigned, by %d\n", BIT_SHIFT);
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr uint16_t byte_shift = BIT_SHIFT >> 3;
        constexpr size_t bytes = sizeof(WORD_TYPE) * len - byte_shift;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // left shift, byte aligned, copy from MSB to LSB

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out) + sizeof(WORD_TYPE) * len;
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in) + bytes;


		u -= sizeof(MachineWord);
        MachineWord x = loadu<MachineWord>(u);
        MachineWord s = op(x);

        // now loop over all remaining full words.  (loop may not execute...)
        for (unsigned int i = 1; i < nMachWord; ++i) {

        	// read next before write - for in-place shift.
        	u -= sizeof(MachineWord);
        	x = loadu<MachineWord>(u);

        	// write
        	w -= sizeof(MachineWord);
            storeu<MachineWord>(w, s);

        	s = op(x);
        }

        // now take care of remainder.
        if (rem > 0) {
        	// read from beginning
        	x = loadu<MachineWord>(in);

        	// write
        	w -= sizeof(MachineWord);
            storeu<MachineWord>(w, s);

        	s = op(x);
        }
        // write out at beginning.
        storeu<MachineWord>(reinterpret_cast<uint8_t*>(out) + byte_shift, s);

        // finally,
        if (byte_shift > 0) memset(out, 0, byte_shift);  // clear the first part

      }
      // right shift byte aligned.  safe for in and out to be same
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
          ((BIT_SHIFT & 0x7) == 0) && (BIT_SHIFT < 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert((-BIT_SHIFT) <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
//    	  printf("right shift transform, large, bytealigned, by %d\n", BIT_SHIFT);
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr uint16_t byte_shift = (-BIT_SHIFT) >> 3;
        constexpr size_t bytes = sizeof(WORD_TYPE) * len - byte_shift;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // shift from LSB to MSB

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in) + byte_shift;

        MachineWord x = loadu<MachineWord>(u);  // faster than memcpy
        MachineWord s = op(x);

        // handle remainders first
        if (rem > 0) {
        	u += rem;

        	// read next before write - for in-place shift.
        	x = loadu<MachineWord>(u);

            storeu<MachineWord>(w, s);
        	w += rem;

        	s = op(x);
        }
    	u += sizeof(MachineWord);


        // now loop over all remaining full words.  (loop may not execute...)
        for (unsigned int i = 1; i < nMachWord; ++i) {
        	// now get current right shifted

        	x = loadu<MachineWord>(u);
        	u += sizeof(MachineWord);

            storeu<MachineWord>(w, s);
            w += sizeof(MachineWord);

            s = op(x);
        }
        storeu<MachineWord>(w, s);

        // finally,
        if (byte_shift > 0) memset(reinterpret_cast<uint8_t *>(out) + bytes, 0, byte_shift);  // clear the last part
      }

      // left shift.  logic slightly different from right shift.  inplace support
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
              ((BIT_SHIFT & 0x7) != 0) && (BIT_SHIFT > 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert(BIT_SHIFT <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        //  	  printf("left shift transform, large,  by %d, bytes %lu, machineword size %lu\n", BIT_SHIFT, len * sizeof(WORD_TYPE), sizeof(MachineWord));
        constexpr uint16_t byte_shift = BIT_SHIFT >> 3;
        constexpr uint16_t bit_shift = BIT_SHIFT & 0x7;
        constexpr size_t byte_len = len * sizeof(WORD_TYPE) - byte_shift;
        constexpr size_t overlap = (sizeof(MachineWord) == 1 ? 0 : 1);
        constexpr size_t nonoverlap = sizeof(MachineWord) - overlap;  //machine word is larger than 1 because of the conditionals


        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = (byte_len - overlap) / nonoverlap;
        constexpr unsigned int rem = (byte_len - overlap) % nonoverlap;

        // left shift, inplace, from MSB to LSB.  need to make sure real data is
        // on the LSB side of the word, so left shift shifts in 0's.  take care of rem last on LSB side.

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out) + len * sizeof(WORD_TYPE) - overlap;
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in) + byte_len - overlap;

		u -= nonoverlap;
        MachineWord x = loadu<MachineWord>(u);
        MachineWord s = bliss::utils::bit_ops::slli<bit_shift>(op(x));

        // now loop over all remaining full words.  (loop may not execute...)
        for (unsigned int i = 1; i < nMachWord; ++i) {

        	// read next before write - for in-place shift.
        	u -= nonoverlap;
        	x = loadu<MachineWord>(u);

        	// write
        	w -= nonoverlap;
            storeu<MachineWord>(w, s);

        	s = bliss::utils::bit_ops::slli<bit_shift>(op(x));
        }

        // now take care of remainder.
        if (rem > 0) {
        	// read from beginning
        	x = loadu<MachineWord>(in);

        	// write
        	w -= nonoverlap;
            storeu<MachineWord>(w, s);

        	s = bliss::utils::bit_ops::slli<bit_shift>(op(x));
        }
        // write out at beginning.
        storeu<MachineWord>(reinterpret_cast<uint8_t*>(out) + byte_shift, s);

        // finally, clear the LSB.
        if (byte_shift > 0) memset(out, 0, byte_shift);
      }

      // right shift. logic slightly different from left shift, inplace support
      template <typename MAX_SIMD_TYPE, int16_t BIT_SHIFT = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
              ((BIT_SHIFT & 0x7) != 0) && (BIT_SHIFT < 0)), int>::type = 1>
      BITS_INLINE void shift_transform(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len], OP const & op) {
    	  static_assert((-BIT_SHIFT) <= (sizeof(WORD_TYPE) << 3), "ERROR: bit shift should not be larger than word size");
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

//  	  printf("right shift transform, large,  by %d, bytes %lu, machineword size %lu\n", BIT_SHIFT, len * sizeof(WORD_TYPE), sizeof(MachineWord));
        constexpr uint16_t byte_shift = (-BIT_SHIFT) >> 3;
        constexpr uint16_t bit_shift = (-BIT_SHIFT) & 0x7;
        constexpr size_t byte_len = len * sizeof(WORD_TYPE) - byte_shift;
        constexpr size_t overlap = (sizeof(MachineWord) == 1 ? 0 : 1);
        constexpr size_t nonoverlap = sizeof(MachineWord) - overlap;  //machine word is larger than 1 because of the conditionals

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = (byte_len - overlap) / nonoverlap;
        constexpr unsigned int rem = (byte_len - overlap) % nonoverlap;

        // right shift, inplace, from LSB to MSB.  need to make sure real data is on the MSB side of the word
        // so right shift shifts in 0s, so take care of rem first on LSB side

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in) + byte_shift;

        MachineWord x = loadu<MachineWord>(u);  // faster than memcpy
        MachineWord s = bliss::utils::bit_ops::srli<bit_shift>(op(x));

        // handle remainders first
        if (rem > 0) {
        	u += rem;

        	// read next before write - for in-place shift.
        	x = loadu<MachineWord>(u);

            storeu<MachineWord>(w, s);
        	w += rem;

        	s = bliss::utils::bit_ops::srli<bit_shift>(op(x));
        }
    	u += nonoverlap;


        // now loop over all remaining full words.  (loop may not execute...)
        for (unsigned int i = 1; i < nMachWord; ++i) {
        	// now get current right shifted

        	x = loadu<MachineWord>(u);
        	u += nonoverlap;

            storeu<MachineWord>(w, s);
            w += nonoverlap;

            s = bliss::utils::bit_ops::srli<bit_shift>(op(x));
        }
        storeu<MachineWord>(w, s);

        if (byte_shift > 0) memset(reinterpret_cast<uint8_t*>(out) + byte_len, 0, byte_shift);  // clear the last part
      }


      //=================== binary comparison ===========================

      // WORD orders bit bitwise transforms should be same between input and output, just shifted
      // positive shift increase value, so shift left.  negative shift shifts right. 0 shift returns original value, not sure if it's no-op.
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE int8_t bit_compare(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

//    	  printf("binary bit compare exact\n");
    	  static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          MachineWord x = loadu<MachineWord>(lhs);
          MachineWord y = loadu<MachineWord>(rhs);
          return (x == y) ? 0 : ((x < y) ? -1 : 1);
      }

      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE int8_t bit_compare(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

    	  static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");

    	  constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        MachineWord x = ::bliss::utils::bit_ops::zero<MachineWord>();
        MachineWord y = ::bliss::utils::bit_ops::zero<MachineWord>();
        load_part<bytes, 0>(x, lhs);
        load_part<bytes, 0>(y, rhs);


//        MachineWord x = 0, y = 0;
//        memcpy(&x, lhs, bytes);
//        memcpy(&y, rhs, bytes);
        return (x == y) ? 0 : ((x < y) ? -1 : 1);
      }

      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)),
              int>::type = 1>
      BITS_INLINE int8_t bit_compare(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

    	  static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // pointers.  from MSB to LSB.
        uint8_t const * u = reinterpret_cast<uint8_t const *>(lhs) + bytes;
        uint8_t const * v = reinterpret_cast<uint8_t const *>(rhs) + bytes;

        MachineWord x, y;

        // no overlap
        for (unsigned int i = 0; i < nMachWord; ++i) {
          // enough bytes.  do an iteration
            u -= sizeof(MachineWord);
            v -= sizeof(MachineWord);

          x = loadu<MachineWord>(u);
          y = loadu<MachineWord>(v);

          if (x != y) return (x < y) ? -1 : 1;
        }

        if (rem > 0) {  // 0 < rem < sizeof(uint64_t)
            x = loadu<MachineWord>(lhs);
            y = loadu<MachineWord>(rhs);

            if (x != y) return (x < y) ? -1 : 1;
        }
        return 0;  // all equal.
      }

      // WORD orders bit bitwise transforms should be same between input and output, just shifted
      // positive shift increase value, so shift left.  negative shift shifts right. 0 shift returns original value, not sure if it's no-op.
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE bool bit_less(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

//        printf("binary bit compare exact\n");
        static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          MachineWord x = loadu<MachineWord>(lhs);
          MachineWord y = loadu<MachineWord>(rhs);
          return (x < y);
      }

      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE bool bit_less(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

        static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        MachineWord x = ::bliss::utils::bit_ops::zero<MachineWord>();
        MachineWord y = ::bliss::utils::bit_ops::zero<MachineWord>();
        load_part<bytes, 0>(x, lhs);
        load_part<bytes, 0>(y, rhs);


//        MachineWord x = 0, y = 0;
//        memcpy(&x, lhs, bytes);
//        memcpy(&y, rhs, bytes);
        return (x < y);
      }

      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)),
              int>::type = 1>
      BITS_INLINE bool bit_less(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

        static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // pointers.  from MSB to LSB.
        uint8_t const * u = reinterpret_cast<uint8_t const *>(lhs) + bytes;
        uint8_t const * v = reinterpret_cast<uint8_t const *>(rhs) + bytes;

        MachineWord x, y;

        // no overlap
        for (unsigned int i = 0; i < nMachWord; ++i) {
          // enough bytes.  do an iteration
            u -= sizeof(MachineWord);
            v -= sizeof(MachineWord);

          x = loadu<MachineWord>(u);
          y = loadu<MachineWord>(v);

          if (x != y) return (x < y);
        }

        if (rem > 0) {  // 0 < rem < sizeof(uint64_t)
            x = loadu<MachineWord>(lhs);
            y = loadu<MachineWord>(rhs);

          return (x < y);
        } else {
          return false;  // all equal.
        }
      }


      // WORD orders bit bitwise transforms should be same between input and output, just shifted
      // positive shift increase value, so shift left.  negative shift shifts right. 0 shift returns original value, not sure if it's no-op.
      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE bool bit_equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

//        printf("binary bit compare exact\n");
        static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          MachineWord x = loadu<MachineWord>(lhs);
          MachineWord y = loadu<MachineWord>(rhs);
          return (x == y);
      }

      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)), int>::type = 1>
      BITS_INLINE bool bit_equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

        static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        MachineWord x = ::bliss::utils::bit_ops::zero<MachineWord>();
        MachineWord y = ::bliss::utils::bit_ops::zero<MachineWord>();
        load_part<bytes, 0>(x, lhs);
        load_part<bytes, 0>(y, rhs);


//        MachineWord x = 0, y = 0;
//        memcpy(&x, lhs, bytes);
//        memcpy(&y, rhs, bytes);
        return (x == y);
      }

      template <typename MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          typename std::enable_if<((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)),
              int>::type = 1>
      BITS_INLINE bool bit_equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {

        static_assert(MAX_SIMD_TYPE::SIMDVal < BIT_REV_SSSE3, "ERROR bit_compare does not support SSSE3 or AVX2");
        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr size_t bytes = sizeof(WORD_TYPE) * len;

        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = bytes / sizeof(MachineWord);
        constexpr unsigned int rem = bytes & (sizeof(MachineWord) - 1);

        // pointers.  from MSB to LSB.
        uint8_t const * u = reinterpret_cast<uint8_t const *>(lhs) + bytes;
        uint8_t const * v = reinterpret_cast<uint8_t const *>(rhs) + bytes;

        MachineWord x, y;

        // no overlap
        for (unsigned int i = 0; i < nMachWord; ++i) {
          // enough bytes.  do an iteration
            u -= sizeof(MachineWord);
            v -= sizeof(MachineWord);

          x = loadu<MachineWord>(u);
          y = loadu<MachineWord>(v);

          if (x != y) return false;
        }

        if (rem > 0) {  // 0 < rem < sizeof(uint64_t)
            x = loadu<MachineWord>(lhs);
            y = loadu<MachineWord>(rhs);

            return (x == y);
        } else {
          return true;  // all equal.
        }
      }


      //========================= reverse + unary transform =================
      /**
       * @brief
       * @details   reverse for a fixed length array.  specialized for the following mutually exclusive cases
       * 				1. len * sizeof(word_type) == sizeof(machineword_type)  1a. shift = 0, 1b shift > 0
       * 				2. len * sizeof(word_type) <  sizeof(machineword_type)  2a. shift = 0, 2b shift > 0
       * 				3. len * sizeof(word_type) >  sizeof(machineword type)  (need to shift bits < sizeof(word_type))
       * 				   note that len * sizeof(word_type) is greater than sizeof(machineword_type) by at least 1 word_type,
       * 				   and since shift is smaller than sizeof(word_type), number of data bits is more than sizeof(machineword type)
       * 				3a.  power of 2 bits and shift%8 == 0  - byte aligned, so no bitwise shift is needed
       * 				3b.  power of 2 bits and shift%8 > 0   - not byte aligned, so bitwise shift is needed.
       * 				3c.  bits = 3 - has an overlap, and shift also needed.  easier to save prev machine word and shift and OR? or re-read from mem?
       *
       * 				NOTE 1a, 1b, 2a, and 2b are explicitly specialized to force compiler NOT to branch.
       * @param out
       * @param in
       * @param op        operator.  example is one that performs reverse and shifting.  another example is reverse and negate.
       * @tparam len      number of words.
       */
      template <typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0, uint16_t BYTE_OVERLAP = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
        		  (PAD_BITS == 0)), uint16_t>::type = 1>  // len is power of 2.
      BITS_INLINE void reverse_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          storeu<MachineWord>(out, op(loadu<MachineWord>(in)));
      }
      template <typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0, uint16_t BYTE_OVERLAP = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) == sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
        		  (PAD_BITS > 0)), uint16_t>::type = 1>  // len is power of 2.
      BITS_INLINE void reverse_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          storeu<MachineWord>(out, bliss::utils::bit_ops::srli<PAD_BITS>(op(loadu<MachineWord>(in))));
      }

      template <typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0, uint16_t BYTE_OVERLAP = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
        		  (PAD_BITS == 0)), uint16_t>::type = 1>  // len is power of 2.
      BITS_INLINE void reverse_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          constexpr size_t bytes = sizeof(WORD_TYPE) * len;

//          MachineWord x = ::bliss::utils::bit_ops::zero<MachineWord>();
//          load_part<bytes, 0>(x, in);
//          x = op(x);

          MachineWord x = op(loadu<MachineWord>(in));
//          memcpy(out, reinterpret_cast<uint8_t *>(&x) + sizeof(MachineWord) - bytes, bytes);
          store_part<bytes, (sizeof(MachineWord) - bytes)>(out, x);
      }

      template <typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0, uint16_t BYTE_OVERLAP = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) < sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
        		  (PAD_BITS > 0)), uint16_t>::type = 1>  // len is power of 2.
      BITS_INLINE void reverse_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {
          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;
          constexpr size_t bytes = sizeof(WORD_TYPE) * len;

          MachineWord x = bliss::utils::bit_ops::srli<PAD_BITS>(op(loadu<MachineWord>(in)));
          store_part<bytes, (sizeof(MachineWord) - bytes)>(out, x);
      }

      template <typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0, uint16_t BYTE_OVERLAP = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
                                           ((PAD_BITS & 0x7) == 0)
                                      ), uint16_t>::type = 1>
      BITS_INLINE void reverse_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {

        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr uint16_t byte_shift = PAD_BITS >> 3;
        constexpr size_t byte_len = len * sizeof(WORD_TYPE) - byte_shift;

        constexpr size_t nonoverlap = sizeof(MachineWord) - BYTE_OVERLAP;

        static_assert(sizeof(MachineWord) > BYTE_OVERLAP, "machineword is BYTE_OVERLAP");
        static_assert(nonoverlap > 0, "overlap is 0 bytes");

        // want to unroll a few iterations (y) (small is better for smaller data), max x (so little waste)
        // so that there is minimal discarded bytes (w) and small overlap bytes (z)
        // (3x + 8w) = ((16-z)y + z) << 3; choose y is 1, z = 1, w = 1.  (i.e, use 15 of 16 bytes)

        // this means that v = (uint8_t*)(out + len) - w = (uint8_t*)(out + len) - 1
        constexpr unsigned int nMachWord = (byte_len - BYTE_OVERLAP) / nonoverlap;
        constexpr size_t rem = (byte_len - BYTE_OVERLAP) % nonoverlap;

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out) + byte_len - BYTE_OVERLAP;
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);

        if (byte_shift > 0) memset(w + BYTE_OVERLAP, 0, byte_shift);  // clear the last part

        for (size_t i = 0; i < nMachWord; ++i) {
          // enough bytes.  do an iteration
          w -= nonoverlap;

          storeu<MachineWord>(w, op(loadu<MachineWord>(u)));

          //          std::cout << std::dec << " shift " << PAD_BITS << " byte_shift " << byte_shift << " bit_shift " << bit_shift << " wordbits " << word_bits << " bytelen " << byte_len << std::endl;
          //          std::cout << "iter " << i << " prev: "; print(tmp);  std::cout << std::endl;
          //          std::cout << "iter " << i << " shift " << std::dec << (word_bits - bit_shift) << " tmp: "; print(tmp);  std::cout << std::endl;
          //          std::cout << "iter " << i << " out: "; print(std::forward<WORD_TYPE const (&)[len]>(out));  std::cout << std::endl;

          u += nonoverlap;
        }

        // take care of remaining bytes.  bytes is >= MachineWord
         if (rem > 0) {   // bytes - 1 + 15
           MachineWord x = op(loadu<MachineWord>(u));
 		   store_part<(rem + BYTE_OVERLAP), sizeof(MachineWord) - (rem + BYTE_OVERLAP)>(out, x);
         }
      }


      template <typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0, uint16_t BYTE_OVERLAP = 0,
          typename WORD_TYPE, size_t len, typename OP,
          typename std::enable_if<(((sizeof(WORD_TYPE) * len) > sizeof(typename MAX_SIMD_TYPE::MachineWord)) &&
                                           ((PAD_BITS & 0x7) > 0)
                                           ), uint16_t>::type = 1>
      BITS_INLINE void reverse_transform(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], OP const & op) {

        using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

        constexpr uint16_t byte_shift = PAD_BITS >> 3;
        constexpr uint16_t bit_shift = PAD_BITS & 0x7;
        constexpr size_t byte_len = len * sizeof(WORD_TYPE) - byte_shift;

        constexpr size_t nonoverlap = sizeof(MachineWord) - BYTE_OVERLAP;
        static_assert(sizeof(MachineWord) >= 3, "machineword is less than 3 bytes");
        static_assert(nonoverlap > 0, "overlap is 0 bytes");

        constexpr uint16_t left_bit_shift = (nonoverlap << 3) - bit_shift;


        // should be at least 1 since shift < sizeof(word), and len * sizeof(word) > sizeof(machword
        constexpr unsigned int nMachWord = (byte_len - BYTE_OVERLAP) / nonoverlap;
        constexpr unsigned int rem = (byte_len - BYTE_OVERLAP) % nonoverlap;

        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out) + byte_len - BYTE_OVERLAP;
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);

        if (byte_shift > 0) memset(w + BYTE_OVERLAP, 0, byte_shift);  // clear the last part

        // no overlap
        MachineWord rev;
        MachineWord prev = bliss::utils::bit_ops::zero<MachineWord>();

        for (unsigned int i = 0; i < nMachWord; ++i) {
          // enough bytes.  do an iteration
          w -= nonoverlap;

          rev = op(loadu<MachineWord>(u));
          storeu<MachineWord>(w,
        		  	  	  	  bliss::utils::bit_ops::bit_or(prev,
                							  	  	  	  	bliss::utils::bit_ops::srli<bit_shift>(rev)));

//          std::cout << std::dec << " shift " << shift << " byte_shift " << byte_shift << " bit_shift " << bit_shift << " wordbits " << word_bits << " bytelen " << byte_len << std::endl;
//          std::cout << "iter " << i << " prev: "; print(tmp);  std::cout << std::endl;
          prev = bliss::utils::bit_ops::slli<left_bit_shift>(rev);  // if bit_shift == 0, tmp will be set to 0.  if WORD - bit_shift == 0, then tmp will be unchanged
//          std::cout << "iter " << i << " shift " << std::dec << (word_bits - bit_shift) << " tmp: "; print(tmp);  std::cout << std::endl;
//          std::cout << "iter " << i << " out: "; print(std::forward<WORD_TYPE const (&)[len]>(out));  std::cout << std::endl;

          u += nonoverlap;
        }

        if (rem > 0) {  // 0 < rem < sizeof(uint64_t)
            rev = op(loadu<MachineWord>(u));
          rev = bliss::utils::bit_ops::bit_or(prev,
                                              bliss::utils::bit_ops::srli<bit_shift>(rev));

           store_part<(rem + BYTE_OVERLAP), sizeof(MachineWord) - (rem + BYTE_OVERLAP)>(out, rev);
        }
      }



     // ================== dynamic sized versions ========================


      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len      number of words.
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<((MAX_SIMD_TYPE == BIT_REV_SWAR) ||
                                           (MAX_SIMD_TYPE == BIT_REV_SEQ)) &&
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len,  uint16_t bit_offset = 0 ) {

        //printf("swar: ");

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= 32, "ERROR: currently reverse does not support 64 BIT_GRUOP_SIZE for SIMD within a register");

        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
        assert(bit_offset == 0);
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) >> 3)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE) << 3) / BIT_GROUP_SIZE) >= 1);

        //memset(out, 0, len);  // needed because we bitwise OR.
        bitgroup_ops<BIT_GROUP_SIZE, MAX_SIMD_TYPE> op64;

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        uint64_t * w = reinterpret_cast<uint64_t *>(out + len);
        uint64_t const * u = reinterpret_cast<uint64_t const *>(in);

        //printf("remainder %lu\n", rem);



        for (; rem >= WordsInUint64; rem -= WordsInUint64) {
          // enough bytes.  do an iteration
          --w;
          *w = op64.reverse(*u);
          ++u;
        }
        if (rem > 0) {  // 0 < rem < sizeof(uint64_t)
          // do another iteration with all the remaining.
          if (len >= WordsInUint64) {  // original length has 8 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            *(reinterpret_cast<uint64_t *>(out)) =
                op64.reverse(*(reinterpret_cast<uint64_t const *>(in + len - WordsInUint64)));
          } else {  // original length is less than 8, so has to do memcpy.  NO CHOICE.
            uint64_t x = op64.reverse(*(reinterpret_cast<uint64_t const *>(in)));
            memcpy(out, reinterpret_cast<WORD_TYPE *>(&x) + (WordsInUint64 - rem), rem * sizeof(WORD_TYPE));
          }
        }

        //printf("done.\n");
        return 0;  // return remainder.
      }



      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<((MAX_SIMD_TYPE == BIT_REV_SWAR) ||
                                           (MAX_SIMD_TYPE == BIT_REV_SEQ)) &&
                                          (BIT_GROUP_SIZE == 3), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len,  uint16_t bit_offset = 0 ) {

        //printf("swar: ");

        size_t bytes = len * sizeof(WORD_TYPE);

        memset(out, 0, bytes);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = bytes;
        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);
        uint16_t init_offset = bit_offset;

        //printf("remainder %lu\n", rem);

        bitgroup_ops<3, MAX_SIMD_TYPE> op64;


        for (; rem >= 8; rem -= 8) {

          // enough bytes.  do an iteration
          w -= 8;
          bit_offset = op64.reverse(w, u, 8, bit_offset);
          u += 8;

          if ((rem > 8) && (bit_offset > 0)) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            w += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.

        }
        if (rem > 0) {
          // do another iteration with all the remaining.
          if (bytes >= 8) {
            // more than 8 in length, so need to compute overlap with previously computed region (bit offsets)
            bit_offset = (((bytes - 8) << 3) - init_offset) % BIT_GROUP_SIZE;
            bit_offset = (bit_offset == 0) ? 0 : (BIT_GROUP_SIZE - bit_offset);
            bit_offset = op64.reverse(reinterpret_cast<uint8_t *>(out),
                                      reinterpret_cast<uint8_t const *>(in) + bytes - 8,
                                      8, bit_offset);

          } else {
            // original length is less than 32, so has to do memcpy
            bit_offset = op64.reverse(reinterpret_cast<uint8_t *>(out),
                                      reinterpret_cast<uint8_t const *>(in) + bytes - rem,
                                      rem, bit_offset);

          }
        }

        return bit_offset;  // return remainder.
      }



//      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE, size_t len>
//      BITS_INLINE typename std::enable_if<(BIT_GROUP_SIZE == 3), uint16_t>::type
//      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], uint16_t bit_offset = 0 ) {
//        switch (bit_offset % 3) {
//          case 0: return reverse<BIT_GROUP_SIZE, MAX_SIMD_TYPE, 0>(::std::forward<WORD_TYPE (&)[len]>(out),
//                                                                   ::std::forward<WORD_TYPE const (&)[len]>(in)); break;
//          case 1: return reverse<BIT_GROUP_SIZE, MAX_SIMD_TYPE, 1>(::std::forward<WORD_TYPE (&)[len]>(out),
//                                                                   ::std::forward<WORD_TYPE const (&)[len]>(in)); break;
//          case 2: return reverse<BIT_GROUP_SIZE, MAX_SIMD_TYPE, 2>(::std::forward<WORD_TYPE (&)[len]>(out),
//                                                                   ::std::forward<WORD_TYPE const (&)[len]>(in)); break;
//          default: return bit_offset;
//            break;
//        }
//      }
//


#ifdef __SSSE3__
      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len      number of words.
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3) &&
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint16_t bit_offset = 0 ) {
        //printf("ssse3: ");
//        if ((sizeof(WORD_TYPE) * len) < sizeof(__m128i)) {
//          return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE *>(out),
//                                                       ::std::forward<WORD_TYPE const *>(in),
//                                                        len,
//                                                        bit_offset);
//        }

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= (sizeof(uint64_t) << 3), "ERROR: currenly reverse does not support 128 BIT_GRUOP_SIZE for SSSE3");

        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
        assert(bit_offset == 0);
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) >> 3 )) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE) << 3) / BIT_GROUP_SIZE) >= 1);

        //memset(out, 0, len);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        __m128i * w = reinterpret_cast<__m128i *>(out + len);
        __m128i const * u = reinterpret_cast<__m128i const *>(in);

        //printf("remainder %lu\n", rem);

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3> op128;


        for (; rem >= WordsInM128; rem -= WordsInM128) {
          // enough bytes.  do an iteration
          --w;
          _mm_storeu_si128(w, op128.reverse(_mm_loadu_si128(u)));
          ++u;
        }
        if (rem > 0) {  // 0 < rem < WordsInM128
          // do another iteration with all the remaining.
          if (len >= WordsInM128) {  // original length has 8 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            _mm_storeu_si128(reinterpret_cast<__m128i *>(out),
                             op128.reverse(_mm_loadu_si128(reinterpret_cast<__m128i const *>(in + len - WordsInM128))));
          } else {  // original length is less than 8, so has to do memcpy.  NO CHOICE.
            __m128i x = op128.reverse(_mm_loadu_si128(reinterpret_cast<__m128i const *>(in)));
            memcpy(out, reinterpret_cast<WORD_TYPE *>(&x) + (WordsInM128 - rem), rem * sizeof(WORD_TYPE));
          }
        }
        //printf("done.\n");

        return 0;  // return remainder.
      }

      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3) &&
                                          (BIT_GROUP_SIZE == 3), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint16_t bit_offset = 0 ) {

        size_t bytes = len * sizeof(WORD_TYPE);


        memset(out, 0, bytes);  // needed because we bitwise OR.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = bytes;
        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);
        uint16_t init_offset = bit_offset;


        bitgroup_ops<3, BIT_REV_SSSE3, false> op128;
        for (; rem >= 16; rem -= 16) {

          // enough bytes.  do an iteration
          w -= 16;
          bit_offset = op128.reverse(w, u, 16, bit_offset);
          u += 16;

          if ((rem > 16) && (bit_offset > 0)) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            w += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
        }
        if (rem > 0) {

          // do another iteration with all the remaining.
          if (bytes >= 16) {
            bit_offset = (((bytes - 16) << 3) - init_offset) % 3;
            bit_offset = (bit_offset == 0) ? 0 : (BIT_GROUP_SIZE - bit_offset);
            bit_offset = op128.reverse(reinterpret_cast<uint8_t *>(out),
                                       reinterpret_cast<uint8_t const *>(in) + bytes - 16,
                                       16, bit_offset);

          } else {
            // original length is less than 32, so has to do memcpy
            bit_offset = op128.reverse(reinterpret_cast<uint8_t *>(out),
                                       reinterpret_cast<uint8_t const *>(in) + bytes - rem,
                                       rem, bit_offset);
          }
        }  // else let SWAR handle it.
        return bit_offset;

      }

#else
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM128 = 16 / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint16_t bit_offset = 0 ) {
        return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE *>(out),
                                                     ::std::forward<WORD_TYPE const *>(in),
                                                      len,
                                                      bit_offset);

      }
#endif


#ifdef __AVX2__
      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len      number of words.
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2) &&
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint16_t bit_offset = 0 ) {

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= (sizeof(__m128i) << 3), "ERROR: currenly reverse does not support 256 BIT_GRUOP_SIZE for SSSE3");



        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
        assert(bit_offset == 0);
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) >> 3)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE) << 3) / BIT_GROUP_SIZE) >= 1);

        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        __m256i * w = reinterpret_cast<__m256i *>(out + len);
        __m256i const * u = reinterpret_cast<__m256i const *>(in);

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2> op256;


        for (; rem >= WordsInM256; rem -= WordsInM256) {
          // enough bytes.  do an iteration
          --w;
          _mm256_storeu_si256(w, op256.reverse(_mm256_loadu_si256(u)));
          ++u;
        }
        if (rem > 0) {  // 0 < rem < WordsInM256
          // do another iteration with all the remaining.
          if (len >= WordsInM256) {  // original length has 8 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out),
                             op256.reverse(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(in + len - WordsInM256))));
          } else {  // original length is less than 8, so has to do memcpy.  NO CHOICE.
            __m256i x = op256.reverse(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(in)));
            memcpy(out, reinterpret_cast<WORD_TYPE *>(&x) + (WordsInM256 - rem), rem * sizeof(WORD_TYPE));
          }
        }

        //printf("done.\n");

        return 0;  // return remainder.
      }

      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2) &&
                                          (BIT_GROUP_SIZE == 3), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint16_t bit_offset = 0 ) {

        size_t bytes = len * sizeof(WORD_TYPE);


        memset(out, 0, bytes);  // needed because we bitwise OR.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = bytes;
        // pointers
        uint8_t * w = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);
        uint16_t init_offset = bit_offset;


        bitgroup_ops<3, BIT_REV_AVX2, false> op256;

        for (; rem >= 32; rem -= 32) {
          // enough bytes.  do an iteration
          w -= 32;
          bit_offset = op256.reverse(w, u, 32, bit_offset);
          u += 32;

          if ((rem > 32) && (bit_offset > 0) ) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            w += 1;
            bit_offset = static_cast<unsigned short>(8) - bit_offset;
          } // else no adjustment is needed.
        }
        if (rem > 0) {
          // do another iteration with all the remaining.   not that offset is tied to the current u position, so can't use memcpy-avoiding approach.
          if (bytes >= 32) {
            bit_offset = (((bytes - 32) << 3) - init_offset) % BIT_GROUP_SIZE;
            bit_offset = (bit_offset == 0) ? 0 : (BIT_GROUP_SIZE - bit_offset);
            bit_offset = op256.reverse(reinterpret_cast<uint8_t *>(out),
                                       reinterpret_cast<uint8_t const *>(in) + bytes - 32,
                                       32, bit_offset);
          } else {
            // original length is less than 32, so has to do memcpy
            bit_offset = op256.reverse(reinterpret_cast<uint8_t *>(out),
                                       reinterpret_cast<uint8_t const *>(in) + bytes - rem,
                                       rem, bit_offset);
          }
        }  // otherwise use either SSSE3 or SWAR
        return bit_offset;

      }


#else
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM256 = 32 / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2), uint16_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint16_t bit_offset = 0 ) {
        // cascade to SSSE3 and let the choice of implementation be decided there.
        return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>( ::std::forward<WORD_TYPE *>(out),
                                                       ::std::forward<WORD_TYPE const *>(in),
                                                        len,
                                                        bit_offset);
      }

#endif





      /**
       * @brief      bitwise negation for a byte array.
       * @details    this is done using 64 bit data types.
       *             if using SIMD, we would need 3 to 4 instructions since there is no native negate
       *             and we would need to load, xor (with FFFF from cmpeq), and store.
       *             this would negate the benefit of using SIMD
       * @param out
       * @param in
       * @param len
       * @return
       */
      template <typename WORD_TYPE>
      BITS_INLINE void
      bit_not(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len) {
        static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be 64 bit or less");

        constexpr unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE);

        // process uint64_t
        size_t len_64 = len - (len & (WordsInUint64 - 1));
        for (size_t i = 0; i < len_64; i += WordsInUint64) {
          *(reinterpret_cast<uint64_t *>(out + i)) = ~(*(reinterpret_cast<uint64_t const *>(in + i)));
        }
        // process remainder.
        for (size_t i = len_64; i < len; ++i) {
          *(out + i) = ~(*(in + i));
        }
      }


      // NEED TO SPECIFY WORD_TYPE and len as template parameters for function when using icc.  auto type deduction has a bug for fixed size arrays..

      //========================== reverse with shift ================

      /**
       * @brief
       * @details   reverse for a fixed length array.  specialized for the following mutually exclusive cases
       * @param out		in and out should NOT be the same array.
       * @param in		in and out should NOT be the same array.
       * @param op        operator.  example is one that performs reverse and shifting.  another example is reverse and negate.
       * @tparam len      number of words.
       * @tparam PAD_BITS  number of 0 padding bits (at MSB end, last word)
       * @return       if the first bitgroup at LSB is a partial group, the number of bits in that partial group.  may be 0.
       */
      template <unsigned int BIT_GROUP_SIZE, typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0,
          typename WORD_TYPE, size_t len>
      BITS_INLINE uint16_t  // len is power of 2.
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

          static_assert(BIT_GROUP_SIZE <= (sizeof(MachineWord) * 4), "ERROR: The maximum bit group size should not be larger than 1/2 of machine word type");
          static_assert(
        		  ((((sizeof(WORD_TYPE) * len << 3) - PAD_BITS) % BIT_GROUP_SIZE) == 0) , "ERROR: bit_group_size needs to divide all bits evenly.");

          ::bliss::utils::bit_ops::bitgroup_ops<BIT_GROUP_SIZE, MAX_SIMD_TYPE::SIMDVal> op;
          // if power of 2 BIT_GROUP_SIZE, then overlap is 0. else there is overlap between successive machineword reverses
          constexpr uint16_t byte_overlap =
              ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0) ? 0 : sizeof(MachineWord) % BIT_GROUP_SIZE;

          bliss::utils::bit_ops::reverse_transform<MAX_SIMD_TYPE, PAD_BITS, byte_overlap, WORD_TYPE, len>(out, in,
                           [&op](MachineWord const & src){ return op.reverse(src); });

          return ((sizeof(WORD_TYPE) * len << 3) - PAD_BITS) % BIT_GROUP_SIZE;
      }


      template <unsigned int BIT_GROUP_SIZE, typename MAX_SIMD_TYPE, uint16_t PAD_BITS = 0,
          typename WORD_TYPE, size_t len>
      BITS_INLINE uint16_t  // len is power of 2.
      reverse_bits_in_byte(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");

          using MachineWord = typename MAX_SIMD_TYPE::MachineWord;

          static_assert(BIT_GROUP_SIZE <= (sizeof(MachineWord) * 4), "ERROR: The maximum bit group size should not be larger than 1/2 of machine word type");
          static_assert(
        		  ((((sizeof(WORD_TYPE) * len << 3) - PAD_BITS) % BIT_GROUP_SIZE) == 0) , "ERROR: bit_group_size needs to divide all bits evenly.");

          ::bliss::utils::bit_ops::bitgroup_ops<BIT_GROUP_SIZE, MAX_SIMD_TYPE::SIMDVal> op;

          bliss::utils::bit_ops::bit_transform<MAX_SIMD_TYPE, WORD_TYPE, len>(out, in,
                           [&op](MachineWord const & src){ return op.reverse_bits_in_byte(src); });

          return ((sizeof(WORD_TYPE) * len << 3) - PAD_BITS) % BIT_GROUP_SIZE;
      }


      //========================== bitwise operations ============
      // no difference between conservative and aggressive.  use specified.

      // bitwise negation
      template <typename MAX_SIMD_TYPE, typename WORD_TYPE, size_t len>
      BITS_INLINE void bit_not(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
//    	  using SIMD_TYPE = BITREV_AUTO_CONSERVATIVE<(len * sizeof(WORD_TYPE)), MAX_SIMD_TYPE>;
    	  using SIMD_TYPE = MAX_SIMD_TYPE;
    	bliss::utils::bit_ops::bit_transform<SIMD_TYPE, WORD_TYPE, len>(out, in,
          [](typename SIMD_TYPE::MachineWord const & src) { return bliss::utils::bit_ops::bit_not(src); });
      }
      // bitwise and
      template <typename MAX_SIMD_TYPE, typename WORD_TYPE, size_t len>
      BITS_INLINE void bit_and(WORD_TYPE (&out)[len], WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
//    	  using SIMD_TYPE = BITREV_AUTO_CONSERVATIVE<(len * sizeof(WORD_TYPE)), MAX_SIMD_TYPE>;
    	  using SIMD_TYPE = MAX_SIMD_TYPE;

    	  bliss::utils::bit_ops::bit_transform<SIMD_TYPE, WORD_TYPE, len>(out, lhs, rhs,
          [](typename SIMD_TYPE::MachineWord const & l,
              typename SIMD_TYPE::MachineWord const & r ) { return bliss::utils::bit_ops::bit_and(l, r); });
      }
      // bitwise or
      template <typename MAX_SIMD_TYPE, typename WORD_TYPE, size_t len>
      BITS_INLINE void bit_or(WORD_TYPE (&out)[len], WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
//    	  using SIMD_TYPE = BITREV_AUTO_CONSERVATIVE<(len * sizeof(WORD_TYPE)), MAX_SIMD_TYPE>;
    	  using SIMD_TYPE = MAX_SIMD_TYPE;

        bliss::utils::bit_ops::bit_transform<SIMD_TYPE, WORD_TYPE, len>(out, lhs, rhs,
          [](typename SIMD_TYPE::MachineWord const & l,
              typename SIMD_TYPE::MachineWord const & r ) { return bliss::utils::bit_ops::bit_or(l, r); });
      }
      // bitwise xor
      template <typename MAX_SIMD_TYPE, typename WORD_TYPE, size_t len>
      BITS_INLINE void bit_xor(WORD_TYPE (&out)[len], WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
//    	  using SIMD_TYPE = BITREV_AUTO_CONSERVATIVE<(len * sizeof(WORD_TYPE)), MAX_SIMD_TYPE>;
    	  using SIMD_TYPE = MAX_SIMD_TYPE;

        bliss::utils::bit_ops::bit_transform<SIMD_TYPE, WORD_TYPE, len>(out, lhs, rhs,
          [](typename SIMD_TYPE::MachineWord const & l,
              typename SIMD_TYPE::MachineWord const & r ) { return bliss::utils::bit_ops::bit_xor(l, r); });
      }


      //========================== shift operations ============
      // memcpy costly.  use conservative to avoid partial memcpy (some times).
      template <typename MAX_SIMD_TYPE, uint16_t BIT_SHIFT, typename WORD_TYPE, size_t len> //,
//       typename ::std::enable_if<(MAX_SIMD_TYPE::SIMDVal <= BIT_REV_SWAR), int>::type = 1 >
      BITS_INLINE void left_shift(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len]) {
//    	  printf("left shifting by %u\n", BIT_SHIFT);
    	  using SIMD_TYPE = BITREV_AUTO_CONSERVATIVE<(len * sizeof(WORD_TYPE)), MAX_SIMD_TYPE>;
    	  bliss::utils::bit_ops::shift_transform<SIMD_TYPE, BIT_SHIFT, WORD_TYPE, len>(out, in,
          [](typename SIMD_TYPE::MachineWord const & src) { return src; });
    	  }
      template <typename MAX_SIMD_TYPE, uint16_t BIT_SHIFT, typename WORD_TYPE, size_t len> //,
//      	  typename ::std::enable_if<(MAX_SIMD_TYPE::SIMDVal <= BIT_REV_SWAR), int>::type = 1 >
      BITS_INLINE void right_shift(WORD_TYPE (&out)[len], WORD_TYPE (&in)[len]) {
//    	  printf("right shifting %lu x %lu byte words by %u\n", len, sizeof(WORD_TYPE), BIT_SHIFT);
    	  using SIMD_TYPE = BITREV_AUTO_CONSERVATIVE<(len * sizeof(WORD_TYPE)), MAX_SIMD_TYPE>;
        bliss::utils::bit_ops::shift_transform<SIMD_TYPE, (0 - (static_cast<int16_t>(BIT_SHIFT))), WORD_TYPE, len>(out, in,
          [](typename SIMD_TYPE::MachineWord const & src) { return src; });
        }

      //========================== comparison operations ===========
      // Aggressive strategy means partial load. but since SWAR, it's okay. (requires memcpy)
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE bool equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
        using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
        return bliss::utils::bit_ops::bit_equal<MAX_SIMD_TYPE, WORD_TYPE, len>(lhs, rhs);
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE bool less(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
        using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
        return bliss::utils::bit_ops::bit_less<MAX_SIMD_TYPE, WORD_TYPE, len>(lhs, rhs);
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE bool greater(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
          using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
          return bliss::utils::bit_ops::bit_less<MAX_SIMD_TYPE, WORD_TYPE, len>(rhs, lhs);
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE bool not_equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
          using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
          return !(bliss::utils::bit_ops::bit_equal<MAX_SIMD_TYPE, WORD_TYPE, len>(lhs, rhs));
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE bool less_equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
          using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
          return !(bliss::utils::bit_ops::bit_less<MAX_SIMD_TYPE, WORD_TYPE, len>(rhs, lhs));
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE bool greater_equal(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
          using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
          return !(bliss::utils::bit_ops::bit_less<MAX_SIMD_TYPE, WORD_TYPE, len>(lhs, rhs));
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE int8_t compare(WORD_TYPE const (&lhs)[len], WORD_TYPE const (&rhs)[len]) {
          using MAX_SIMD_TYPE = BITREV_AUTO_AGGRESSIVE<(len * sizeof(WORD_TYPE)), BITREV_SWAR>;
          return bliss::utils::bit_ops::bit_compare<MAX_SIMD_TYPE, WORD_TYPE, len>(lhs, rhs);
      }
    } // namespace bit_ops


  } // namespace utils

} // namespace bliss



#if defined __GNUC__ && __GNUC__>=6
// disable __m128i and __m256i ignored attribute warning in gcc
  #pragma GCC diagnostic pop
#endif





#endif /* SRC_UTILS_BITGROUP_OPS_HPP_ */
