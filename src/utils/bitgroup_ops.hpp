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

#include "utils/logging.h"
#include "bliss-config.hpp"
#ifdef USE_SIMD
#include <x86intrin.h>   // all intrinsics.  will be enabled based on compiler flag such as __SSSE3__ internally.
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

#define _bliss_bswap_16(x)  ((x >>  8) &  0xFF) | ((x &  0xFF) <<  8);

// TODO: reverse_groups vs reverse_in_group


namespace bliss {

  namespace utils {

    namespace bit_ops {

      static constexpr unsigned char BIT_REV_SEQ = 0;
      static constexpr unsigned char BIT_REV_SWAR = 1;   // SIMD Within A Register
      static constexpr unsigned char BIT_REV_SSSE3 = 2;
      static constexpr unsigned char BIT_REV_AVX2 = 4;




      /// shift right by number of bits less than 8.
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE srli(WORD_TYPE const val, uint_fast8_t shift) {
        return val >> shift;
      }
      /// shift right by number of bits less than 8.
      template <typename WORD_TYPE, uint8_t shift>
      BITS_INLINE WORD_TYPE srli(WORD_TYPE const val) {
        return val >> shift;
      }
      /// shift left by number of bits less than 8.
      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE slli(WORD_TYPE const val, uint_fast8_t shift) {
        return val << shift;
      }
      /// shift left by number of bits less than 8.
      template <typename WORD_TYPE, uint8_t shift>
      BITS_INLINE WORD_TYPE slli(WORD_TYPE const val) {
        return val << shift;
      }

      template <typename WORD_TYPE>
      BITS_INLINE WORD_TYPE negate(WORD_TYPE const & u) {
        return ~u;
      }

      /**
       * @brief base bit reverse type. base template only
       * @tparam WORD_TYPE          type of input word
       * @tparam BIT_GROUP_SIZE     number of bits in a group to be reversed.  supports any value less than 8, and powers of 2.  tested 3 and powers of 2.  cannot exceed word size.
       * @tparam BIT_REV_SIMD_TYPE  type of algorithm to use based on available hardware
       * @tparam DUMMY              here so that we can define static constexpr array variables in the header file.
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char BIT_REV_SIMD_TYPE = BIT_REV_SEQ, bool POW2 = ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0)>
      struct bitgroup_ops {
          //  default imple is for SEQuential bit reverse.  this is defined for all bit_group_sizes.
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < (sizeof(uint64_t) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in uint64_t");
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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
            if ((len * 8) < BIT_GROUP_SIZE) return bit_offset; //throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            // if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be less than 8.  for larger, shift input instead.");
            assert(bit_offset < 8);

//            if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//              throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
            assert((len % ((BIT_GROUP_SIZE + 7) / 8)) == 0);

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

            return (len * 8 - bit_offset) % BIT_GROUP_SIZE;
          }

          /**
           * @brief  reverse the lower portion of WORD, starting from 0, up to the largest multiples of BIT_GROUP_SIZE that still fits.
           *
           * @param v    input word.  LSB aligned to start of a bit group.
           * @return     bit reversed word.  MSB aligned to end of the originally lowest bit group.
           */
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS >= (sizeof(WORD_TYPE) * 8)), int>::type = 1>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & v, uint8_t bit_offset = 0) {
//            static_assert(BIT_GROUP_SIZE == (sizeof(WORD_TYPE) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");
            return v;
          }

          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS < (sizeof(WORD_TYPE) * 8)), int>::type = 1>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & v, uint8_t bit_offset = 0) {

            // copy the source data
            WORD_TYPE w = 0;
            WORD_TYPE u = v;

            unsigned int rem = (sizeof(WORD_TYPE) * 8 - bit_offset) % BIT_GROUP_SIZE;
            //WORD_TYPE rem_mask = static_cast<WORD_TYPE>(~(std::numeric_limits<WORD_TYPE>::max() >> rem));

            u >>= bit_offset;

            // now do the middle part
            for (size_t i = bit_offset; i < (sizeof(WORD_TYPE) * 8 - rem); i += BIT_GROUP_SIZE) {
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


      };

      /// partial template specialization for SWAR based bit reverse.  uses BSWAP when appropriate, else use bit shift  this is defined only for bit_group_sizes that are powers of 2 up to sizeof WORD_TYPE)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SWAR, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < (sizeof(uint64_t) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in uint64_t");
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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
            if ( (len * 8) < BIT_GROUP_SIZE ) return bit_offset; // throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be 0, since we are ORing the first and last bytes.
            //if (bit_offset > 0) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be 0 for SWAR with power of 2 BIT_GROUP_SIZE.");
            assert(bit_offset == 0 );

            //if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
            //  throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
            assert((len % ((BIT_GROUP_SIZE + 7) / 8)) == 0);


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
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS > (sizeof(WORD_TYPE) * 8)), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) {
            return u;
          }
          // bitgroup size is same size as the word type in bits.
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS == (sizeof(WORD_TYPE) * 8)), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) {
            return u;
          }
          // this specialization will work for all WORD_TYPE sizes.
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS <= 8), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) {
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
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS < (sizeof(WORD_TYPE) * 8)) && (sizeof(WORD_TYPE) == 8), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) {
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
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS > 8) && (BITS < (sizeof(WORD_TYPE) * 8)) && (sizeof(WORD_TYPE) == 4), WORD_TYPE>::type
          reverse(WORD_TYPE const &u) {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            return ((u >> 16) & mask16) | ((u & mask16) << 16);
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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
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
                  v7i &= (mask_all >> (64 - len * 8));
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


            return (len * 8 - bit_offset) % 3;
          }

          // dispatcher
          template <typename WORD_TYPE>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const & u, uint8_t const & bit_offset) {
            switch (bit_offset % 3) {
              case 0: return reverse<WORD_TYPE, 0>(::std::forward<WORD_TYPE const>(u)); break;
              case 1: return reverse<WORD_TYPE, 1>(::std::forward<WORD_TYPE const>(u)); break;
              case 2: return reverse<WORD_TYPE, 2>(::std::forward<WORD_TYPE const>(u)); break;
              default: return 0; break;
            }
          }

          /// reverse in groups of 3 bits.  NOTE: within the offset or remainder, the middle bit remains,
          /// while the high and low bits are zeroed during the reversal.  this makes OR'ing with adjacent
          /// entries simple.
          template <typename WORD_TYPE, uint8_t offset = 0>
          BITS_INLINE WORD_TYPE reverse(WORD_TYPE const &u) {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            // and swap the 1st and 3rd elements in the group.
            // we need 3 patterns for shifting based on the offset.
            // offset is at the LSB position of the ORIGINAL word

            WORD_TYPE v = u;
            // 2 shifts, 2 ors, 3 ands
            //switch ((sizeof(WORD_TYPE) * 8 - offset) % 3) {
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
      BITS_INLINE __m128i srli(__m128i val, uint_fast8_t shift) {
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

        if (shift > 64)
        	return _mm_srli_epi64(tmp, shift - 64);

        // 64 bit exactly
        return tmp;
      }

      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
      /// DO THIS BECAUSE BITWISE SHIFT does not cross the epi64 boundary.
      template <uint8_t shift>
      BITS_INLINE __m128i srli(__m128i val) {
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

        if (shift > 64)
          return _mm_srli_epi64(tmp, shift - 64);

        // 64 bit exactly
        return tmp;
      }

      /// shift left by number of bits less than 8.
      template <>
      BITS_INLINE __m128i slli(__m128i val, uint_fast8_t shift) {
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
      template <uint8_t shift>
      BITS_INLINE __m128i slli(__m128i val) {
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


      template <>
      BITS_INLINE __m128i negate(__m128i const & u) {
        return _mm_xor_si128(u, _mm_cmpeq_epi8(u, u));  // no native negation operator, so use xor
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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
            if ((len * 8) < BIT_GROUP_SIZE) return bit_offset; // throw ::std::invalid_argument("ERROR reversing byte array: length is 0");
            //if (len > 16) throw std::invalid_argument("ERROR:  length of array should be <= 16.  For longer, use a different bitgroup_ops (AVX2) or break it up.");
            assert ( len <= 16);
//            if (((len * 8) % BIT_GROUP_SIZE) > 0)
//              throw ::std::invalid_argument("ERROR reversing byte array:  len in bits needs to be a multiple of BIT_GROUP_SIZE");
            assert(((len * 8) % BIT_GROUP_SIZE) == 0);

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
          reverse(__m128i const & u) {
            static_assert(BIT_GROUP_SIZE == 128, "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }

          // groupsize =1 version:  4 loads, 3 shuffles, 3 logical ops, 1 shift op
          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS <= 8), __m128i>::type
          reverse(__m128i const & u) {

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
          reverse(__m128i const & u) {

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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
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
            return (len * 8 - bit_offset) % 3;

          }


          BITS_INLINE __m128i reverse(__m128i const & u, uint8_t bit_offset = 0 ) {
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
            template <uint8_t offset = 0>
            BITS_INLINE __m128i reverse(__m128i const & u ) {

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

      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
      template <>
      BITS_INLINE __m256i srli(__m256i val, uint_fast8_t shift) {
    	  if (shift == 0) return val;
    	  if (shift >= 256) return _mm256_setzero_si256();

    	// goal:  get next higher 64 bits to current position.
    	  // conditionals:  shift byt 0 to 63, 64 to 128, 128 to 192, and 192 to 256
    	  // code below can be modified to shift up to 64 bits, and with conditional, up to 128 bits?
    	  // 64, 192: p+alignr is enough
    	  // 128: permute is enough
    	  // 0-63:  same as below - 5 instructions
    	  // 65-127, 129-191: permute + 2 alignr, then or(slli, srli)  total 6 instructions
    	  // 193-255: permute + alignr, then right shift. 3 instructions

    	  //alternatives:  permute4x64 (all 256bits, but can't zero), s(l/r)li_si256 (128 bit lanes)
    	 // alternative, do the 64 bit shift, then do remaining.

        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
        //                        to shift right, b = val, N = 1, lower lane of a should be upper lane of val, upper lane of a should be 0.  possibly 2 ops. permute latency is 3
        __m256i hi = _mm256_permute2x128_si256(val, val, 0x83);  // 1000.0011 higher lane is 0, lower lane is higher lane of val
        if (shift == 128) return hi;

        uint_fast8_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.
        __m256i v1 = _mm256_alignr_epi8(hi, val, byte_shift); // alignr:  src1 low OR src2 low, then right shift bytes, put into output low.  same with high
        uint_fast8_t bit64_shift = shift & 0x3F;        // shift within the 64 bit block

        if (bit64_shift == 0) return v1;  // exact byte alignment.  return it.
        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
        //    0   hi
        //    hi  lo    => dest_lo =  hi lo shift
        //                 dest_hi =   0 hi shift
        if (byte_shift == 24) return _mm256_srli_epi64(v1, bit64_shift);

        // else we need one more alignr
        __m256i v2 = _mm256_alignr_epi8(hi, val, byte_shift + 8);

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        return _mm256_or_si256(_mm256_slli_epi64(v2, 64 - bit64_shift), _mm256_srli_epi64(v1, bit64_shift));
      }

      /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
      template <uint8_t shift>
      BITS_INLINE __m256i srli(__m256i val) {
    	  if (shift == 0) return val;
    	  if (shift >= 256) return _mm256_setzero_si256();

    	// goal:  get next higher 64 bits to current position.
    	  // conditionals:  shift byt 0 to 63, 64 to 128, 128 to 192, and 192 to 256
    	  // code below can be modified to shift up to 64 bits, and with conditional, up to 128 bits?
    	  // 64, 192: p+alignr is enough
    	  // 128: permute is enough
    	  // 0-63:  same as below - 5 instructions
    	  // 65-127, 129-191: permute + 2 alignr, then or(slli, srli)  total 6 instructions
    	  // 193-255: permute + alignr, then right shift. 3 instructions

    	  //alternatives:  permute4x64 (all 256bits, but can't zero), s(l/r)li_si256 (128 bit lanes)
    	 // alternative, do the 64 bit shift, then do remaining.

        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
        //                        to shift right, b = val, N = 1, lower lane of a should be upper lane of val, upper lane of a should be 0.  possibly 2 ops. permute latency is 3
        __m256i hi = _mm256_permute2x128_si256(val, val, 0x83);  // 1000.0011 higher lane is 0, lower lane is higher lane of val
        if (shift == 128) return hi;

        constexpr uint_fast8_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.
        __m256i v1 = _mm256_alignr_epi8(hi, val, byte_shift); // alignr:  src1 low OR src2 low, then right shift bytes, put into output low.  same with high

        constexpr uint_fast8_t bit64_shift = shift & 0x3F;        // shift within the 64 bit block

        if (bit64_shift == 0) return v1;  // exact byte alignment.  return it.
        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
        //    0   hi
        //    hi  lo    => dest_lo =  hi lo shift
        //                 dest_hi =   0 hi shift
        if (byte_shift == 24) return _mm256_srli_epi64(v1, bit64_shift);

        // else we need one more alignr
        __m256i v2 = _mm256_alignr_epi8(hi, val, byte_shift + 8);

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        return _mm256_or_si256(_mm256_slli_epi64(v2, 64 - bit64_shift), _mm256_srli_epi64(v1, bit64_shift));
      }



      /// shift left by number of bits less than 8.
      template <>
      BITS_INLINE __m256i slli(__m256i val, uint_fast8_t shift) {
    	  if (shift == 0) return val;
    	  if (shift >= 256) return _mm256_setzero_si256();

    	  // val is  A B C D

        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
        //                        to shift left, a = val, N = 15, lower lane of b should be 0, higher lane of b should be lower lane of val.  possibly 2 ops, permute latency is 3.
        __m256i lo = _mm256_permute2x128_si256(val, val, 0x08);  // lower lane is 0, higher lane is lower lane of val
        if (shift == 128) return lo;   // lo is C D 0 0

        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
        //    hi  lo
        //    lo   0    => dest_lo =  hi lo shift
        //                 dest_hi =   0 hi shift
        //uint_fast8_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.

        // TODO - left shift is trickier.
		__m256i v1 = (shift < 128) ? _mm256_alignr_epi8(val, lo, 8) :  // permute essentially shifted by 16
				//  else shift > 128
									 _mm256_permute4x64_epi64(lo, 0x80); // permute and shift lo left more.  2 0 0 0 -> 0x80
		// if shift < 128, v1 is B C D 0.  else v1 is D 0 0 0

        uint_fast8_t bit64_shift = shift & 0x3F;        // shift within the 64 bit block

        if (bit64_shift == 0) return v1;

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        __m256i lv, rv;
        switch (shift >> 6) {
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
        return _mm256_or_si256(_mm256_srli_epi64(lv, 64 - shift), _mm256_slli_epi64(rv, bit64_shift));
      }


      /// shift left by number of bits less than 8.
      template <uint8_t shift>
      BITS_INLINE __m256i slli(__m256i val) {
    	  if (shift == 0) return val;
    	  if (shift >= 256) return _mm256_setzero_si256();

    	  // val is  A B C D

        // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
        // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
        //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
        //                        to shift left, a = val, N = 15, lower lane of b should be 0, higher lane of b should be lower lane of val.  possibly 2 ops, permute latency is 3.
        __m256i lo = _mm256_permute2x128_si256(val, val, 0x08);  // lower lane is 0, higher lane is lower lane of val
        if (shift == 128) return lo;   // lo is C D 0 0

        // together, permute puts the high lane in low of tmp, then tmp and val's low parts are shifted, and high parts are shifted.
        //    hi  lo
        //    lo   0    => dest_lo =  hi lo shift
        //                 dest_hi =   0 hi shift
        //uint_fast8_t byte_shift = (shift & 0xC0) >> 3;  // multiple of 8 bytes, so that we can avoid some ors.

        // TODO - left shift is trickier.
		__m256i v1 = (shift < 128) ? _mm256_alignr_epi8(val, lo, 8) :  // permute essentially shifted by 16
				//  else shift > 128
									 _mm256_permute4x64_epi64(lo, 0x80); // permute and shift lo left more.  2 0 0 0 -> 0x80
		// if shift < 128, v1 is B C D 0.  else v1 is D 0 0 0

        constexpr uint_fast8_t bit64_shift = shift & 0x3F;        // shift within the 64 bit block

        if (bit64_shift == 0) return v1;

        // then left shift to get the bits in the right place
        // shift the input value via epi64 by the number of shifts
        // then or together.
        __m256i lv, rv;
        switch (shift >> 6) {
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
        return _mm256_or_si256(_mm256_srli_epi64(lv, 64 - shift), _mm256_slli_epi64(rv, bit64_shift));
      }

      template <>
      BITS_INLINE __m256i negate(__m256i const & u) {
        return _mm256_xor_si256(u, _mm256_cmpeq_epi8(u, u));  // no native negation operator
      }

      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 256bit)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE < 256, "ERROR: BIT_GROUP_SIZE is greater than number of bits in __m256i");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");

          const __m256i rev_idx_lane8;
          const __m256i rev_idx_lane16;
          const __m256i rev_idx32;
          const __m256i _mask_lo;
          const __m256i lut1_lo;
          const __m256i lut1_hi;
          const __m256i lut2_lo;
          const __m256i lut2_hi;
          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_AVX2;

          bitgroup_ops() :
            rev_idx_lane8( _mm256_setr_epi32(0x0C0D0E0F, 0x08090A0B, 0x04050607, 0x00010203, 0x0C0D0E0F, 0x08090A0B, 0x04050607, 0x00010203)),
            rev_idx_lane16(_mm256_setr_epi32(0x0D0C0F0E, 0x09080B0A, 0x05040706, 0x01000302, 0x0D0C0F0E, 0x09080B0A, 0x05040706, 0x01000302)),
            rev_idx32(     _mm256_setr_epi32(0x00000007, 0x00000006, 0x00000005, 0x00000004, 0x00000003, 0x00000002, 0x00000001, 0x00000000)),
            _mask_lo(      _mm256_setr_epi32(0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F)),
            lut1_lo(       _mm256_setr_epi32(0x0c040800, 0x0e060a02, 0x0d050901, 0x0f070b03, 0x0c040800, 0x0e060a02, 0x0d050901, 0x0f070b03)),
            lut1_hi(       _mm256_setr_epi32(0xc0408000, 0xe060a020, 0xd0509010, 0xf070b030, 0xc0408000, 0xe060a020, 0xd0509010, 0xf070b030)),
            lut2_lo(       _mm256_setr_epi32(0x0c080400, 0x0d090501, 0x0e0a0602, 0x0f0b0703, 0x0c080400, 0x0d090501, 0x0e0a0602, 0x0f0b0703)),
            lut2_hi(       _mm256_setr_epi32(0xc0804000, 0xd0905010, 0xe0a06020, 0xf0b07030, 0xc0804000, 0xd0905010, 0xe0a06020, 0xf0b07030))
            {
          }

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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
            if ((len * 8) < BIT_GROUP_SIZE) return bit_offset; // throw ::std::invalid_argument("ERROR reversing byte array: length is 0");
            //if (len > 32) throw std::invalid_argument("ERROR:  length of array should be <= 32.  For longer, break it up.");
            assert(len <= 32);

//            if (((len * 8) % BIT_GROUP_SIZE) > 0)
//              throw ::std::invalid_argument("ERROR reversing byte array:  len in bits needs to be a multiple of BIT_GROUP_SIZE");
            assert(((len * 8) % BIT_GROUP_SIZE) == 0);

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
          reverse(__m256i const & u) {
            static_assert(BIT_GROUP_SIZE == 256, "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE>
          BITS_INLINE typename std::enable_if<(BITS <= 8), __m256i>::type
          reverse(__m256i const & u) {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            __m256i v = u;

            // if bit group is <= 8, then need to shuffle bytes, and then shuffle bits.
            v = _mm256_shuffle_epi8(v, rev_idx_lane8);     // reverse 8 in each of 128 bit lane          // AVX2
            v = _mm256_permute2x128_si256(v, v, 0x03);							    // then swap the lanes                           // AVX2.  latency = 3

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
          reverse(__m256i const & u) {


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

      };


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

          const __m256i mask3lo;
          const __m256i mask3mid;
          const __m256i mask3hi;


          bitgroup_ops() :
            mask3lo( _mm256_setr_epi32(0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492)),
            mask3mid(_mm256_setr_epi32(0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924)),
            mask3hi( _mm256_setr_epi32(0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249, 0x92492492, 0x24924924, 0x49249249))
          {}


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
          BITS_INLINE uint8_t reverse(uint8_t * out, uint8_t const * in, size_t const & len, uint8_t const bit_offset = 0) {
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
            return (len * 8 - bit_offset) % 3;
          }

          BITS_INLINE __m256i reverse(__m256i const & u, uint8_t bit_offset = 0 ) {
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
          template <uint8_t offset = 0>
          BITS_INLINE __m256i reverse(__m256i const & u) {

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
              v = _mm256_or_si256(srli(mid, 2), _mm256_or_si256(lo, slli(hi, 2)) );
              break;
            case 1:
              // rem == 1:                    hi, lo, mid
                v = _mm256_or_si256(srli(lo, 2), _mm256_or_si256(hi, slli(mid, 2)) );
              break;
            default:
              // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
                v = _mm256_or_si256(srli(hi, 2), _mm256_or_si256(mid, slli(lo, 2)) );
              break;
            }
            // alternatively,  for cross-boundary shift, do byte shift, and reverse shift (8-bit_offset), then or. with shifted (offset)->  each CBS is finished in 4 ops w/o load.

            //========================== first reverse bits in groups of 1.
            return bit_rev_1.reverse(v);
            //=========================== done reverse bits in groups of 1.

          }


      };

      template <bool POW2>
      constexpr uint8_t bitgroup_ops<3, BIT_REV_AVX2, POW2>::BLISS_ALIGNED_ARRAY(mask_all, 64, 32);



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
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len,  uint8_t bit_offset = 0 ) {

        //printf("swar: ");

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= 32, "ERROR: currently reverse does not support 64 BIT_GRUOP_SIZE for SIMD within a register");

        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
        assert(bit_offset == 0);
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) == 0);

        //memset(out, 0, len);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        uint64_t * v = reinterpret_cast<uint64_t *>(out + len);
        uint64_t const * u = reinterpret_cast<uint64_t const *>(in);

        //printf("remainder %lu\n", rem);

        bitgroup_ops<BIT_GROUP_SIZE, MAX_SIMD_TYPE> op64;


        for (; rem >= WordsInUint64; rem -= WordsInUint64) {
          // enough bytes.  do an iteration
          --v;
          *v = op64.reverse(*u);
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

      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len      number of words.
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE, size_t len,
          unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<((MAX_SIMD_TYPE == BIT_REV_SWAR) ||
                                           (MAX_SIMD_TYPE == BIT_REV_SEQ)) &&
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {

        //printf("swar: ");

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= 32, "ERROR: currently reverse does not support 64 BIT_GRUOP_SIZE for SIMD within a register");

        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        static_assert(((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) == 0,
                      "ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");

        //memset(out, 0, len);  // needed because we bitwise OR.

        // pointers
        uint64_t * v = reinterpret_cast<uint64_t *>(out + len);
        uint64_t const * u = reinterpret_cast<uint64_t const *>(in);

        //printf("remainder %lu\n", rem);

        bitgroup_ops<BIT_GROUP_SIZE, MAX_SIMD_TYPE> op64;

        for (size_t i = 0; i < (len / WordsInUint64); ++i) {
          // enough bytes.  do an iteration
          --v;
          *v = op64.reverse(*u);
          ++u;
        }

        constexpr size_t rem = (len % WordsInUint64);
        if (rem > 0) {  // 0 < rem < sizeof(uint64_t)
          // do another iteration with all the remaining.
          if (len > WordsInUint64) {  // original length has 8 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
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
                                          (BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len,  uint8_t bit_offset = 0 ) {

        //printf("swar: ");

        size_t bytes = len * sizeof(WORD_TYPE);

        memset(out, 0, bytes);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = bytes;
        // pointers
        uint8_t * v = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);
        uint8_t init_offset = bit_offset;

        //printf("remainder %lu\n", rem);

        bitgroup_ops<3, MAX_SIMD_TYPE> op64;


        for (; rem >= 8; rem -= 8) {

          // enough bytes.  do an iteration
          v -= 8;
          bit_offset = op64.reverse(v, u, 8, bit_offset);
          u += 8;

          if ((rem > 8) && (bit_offset > 0)) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.

        }
        if (rem > 0) {
          // do another iteration with all the remaining.
          if (bytes >= 8) {
            // more than 8 in length, so need to compute overlap with previously computed region (bit offsets)
            bit_offset = ((bytes - 8) * 8 - init_offset) % BIT_GROUP_SIZE;
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

      // SWAR reverse for fixed size arrays, when bit group size is 3.
      // not that offset is always going to be 0.
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<((MAX_SIMD_TYPE == BIT_REV_SWAR) ||
                                           (MAX_SIMD_TYPE == BIT_REV_SEQ)) &&
                                          (BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {

        // because we have a fixed size array, we can compute the offsets for each WORD exactly.
        // we can then directly call the right version of
        constexpr size_t bytes = len * sizeof(WORD_TYPE);

        //memset(out, 0, bytes);  // needed because we bitwise OR.
        // decide which one to use.

        // want to unroll a few iterations (y) (small is better for smaller data),
        // so that there is minimal discarded bytes (w) and small overlap bytes (z)
        // (3x + 8w) = ((8-z)y + z) * 8; smallest y is 2, z = 1, w = 0.
        // example of z = w = 2 - only 6 of 8 bytes are useful, 25% waste.  so pick y = 2 case.

        // this means that v = (uint8_t*)(out + len) - w = (uint8_t*)(out + len)


        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        // pointers
        uint8_t * v = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);

        bitgroup_ops<3, MAX_SIMD_TYPE> op64;

        // TODO: any reason to do it as a 1 bit reversal followed by an in-group reverse?

        uint64_t tmp;

        // if 64 bit, then 1 bit remains.  2 iterations consumes 15 bytes, which is a multiple of 3,
        // so offset is back to 0.  we can use this fact to partially unroll loops.
        for (size_t i = 0; i < (bytes / 15); i += 15) {
          // enough bytes.  do an iteration
          v -= 15;
          *(reinterpret_cast<uint64_t *>(v)) =
              op64.reverse<1>(*(reinterpret_cast<uint64_t const *>(u + 7)));
          tmp = static_cast<uint64_t>(*(v + 7));
          *(reinterpret_cast<uint64_t *>(v + 7)) =
              tmp | op64.reverse<0>(*(reinterpret_cast<uint64_t const *>(u)));
          u += 15;
        }
        // remaining bytes:
        if ((bytes % 15) >= 8) {
          *(reinterpret_cast<uint64_t *>(v - 8)) =
              op64.reverse<0>(*(reinterpret_cast<uint64_t const *>(u)));
        }

        constexpr size_t rem = (bytes % 15) % 8;
        if (rem > 0) {
          if (bytes >= sizeof(uint64_t)) {  // there are enough bytes. we have at least 1 bytes that completely overlaps
            // we can just keep 1 completed byte from previous reverse call.
            tmp = *(reinterpret_cast<uint64_t *>(out)) & 0xFF00000000000000; // get an internal byte from prev reverse.
            *(reinterpret_cast<uint64_t *>(out)) = tmp |
                op64.reverse<(2 - ((bytes - 7) * 8) % 3)>(*(reinterpret_cast<uint64_t const *>(in + len) - 1));
          } else { // too short
            tmp = op64.reverse<0>(*(reinterpret_cast<uint64_t const *>(in)));
            memcpy(out, reinterpret_cast<uint8_t *>(&tmp) + (sizeof(uint64_t) - bytes), bytes);
          }
        }

        return (len * sizeof(WORD_TYPE) * 8) % 3;  // return remainder.
      }


//      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE, size_t len>
//      BITS_INLINE typename std::enable_if<(BIT_GROUP_SIZE == 3), uint8_t>::type
//      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len], uint8_t bit_offset = 0 ) {
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
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint8_t bit_offset = 0 ) {
        //printf("ssse3: ");
        if ((sizeof(WORD_TYPE) * len) <= sizeof(__m128i)) {
          return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE *>(out),
                                                       ::std::forward<WORD_TYPE const *>(in),
                                                        len,
                                                        bit_offset);
        }

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= sizeof(uint64_t) * 8, "ERROR: currenly reverse does not support 128 BIT_GRUOP_SIZE for SSSE3");

        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
        assert(bit_offset == 0);
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) == 0);

        //memset(out, 0, len);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        __m128i * v = reinterpret_cast<__m128i *>(out + len);
        __m128i const * u = reinterpret_cast<__m128i const *>(in);

        //printf("remainder %lu\n", rem);

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3> op128;


        for (; rem >= WordsInM128; rem -= WordsInM128) {
          // enough bytes.  do an iteration
          --v;
          _mm_storeu_si128(v, op128.reverse(_mm_loadu_si128(u)));
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
      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len      number of words.
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE,
        typename WORD_TYPE, size_t len,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3) &&
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
        //printf("ssse3: ");
        if ((sizeof(WORD_TYPE) * len) < sizeof(__m128i)) {
          return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE (&)[len]>(out),
                                                       ::std::forward<WORD_TYPE const (&)[len]>(in));
        }

        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= sizeof(uint64_t) * 8, "ERROR: currenly reverse does not support 128 BIT_GRUOP_SIZE for SSSE3");

        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) == 0);

        //memset(out, 0, len);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        // pointers
        __m128i * v = reinterpret_cast<__m128i *>(out + len);
        __m128i const * u = reinterpret_cast<__m128i const *>(in);

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3> op128;


        for (size_t i = 0; i < (len / WordsInM128);  ++i) {
          // enough bytes.  do an iteration
          --v;
          _mm_storeu_si128(v, op128.reverse(_mm_loadu_si128(u)));
          ++u;
        }

        constexpr size_t rem = (len % WordsInM128);
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
        return 0;  // return remainder.
      }
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3) &&
                                          (BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint8_t bit_offset = 0 ) {

//        if ((sizeof(WORD_TYPE) * len) < sizeof(__m128i)) {
//          // SSSE3 right now is most performant for groups size of 3.  so always use it.
//          return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>( ::std::forward<WORD_TYPE *>(out),
//                                                         ::std::forward<WORD_TYPE const *>(in),
//                                                          len,
//                                                          bit_offset);
//        }
        size_t bytes = len * sizeof(WORD_TYPE);


        memset(out, 0, bytes);  // needed because we bitwise OR.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = bytes;
        // pointers
        uint8_t * v = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);
        uint8_t init_offset = bit_offset;


        bitgroup_ops<3, BIT_REV_SSSE3, false> op128;
        for (; rem >= 16; rem -= 16) {

          // enough bytes.  do an iteration
          v -= 16;
          bit_offset = op128.reverse(v, u, 16, bit_offset);
          u += 16;

          if ((rem > 16) && (bit_offset > 0)) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
        }
        if (rem > 0) {

          // do another iteration with all the remaining.
          if (bytes >= 16) {
            bit_offset = ((bytes - 16) * 8 - init_offset) % 3;
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

      // SSSE3 reverse for fixed size arrays, when bit group size is 3.
      // not that offset is always going to be 0.
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3) &&
                                          (BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {

        if ((sizeof(WORD_TYPE) * len) < sizeof(__m128i)) {
          // SSSE3 right now is most performant for groups size of 3.  so always use it.
          return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE (&)[len]>(out),
                                                       ::std::forward<WORD_TYPE const (&)[len]>(in));
        }
        constexpr size_t bytes = len * sizeof(WORD_TYPE);


        //memset(out, 0, bytes);  // needed because we bitwise OR.

        // want to unroll a few iterations (y) (small is better for smaller data), max x (so little waste)
        // so that there is minimal discarded bytes (w) and small overlap bytes (z)
        // (3x + 8w) = ((16-z)y + z) * 8; choose y is 1, z = 1, w = 1.  (i.e, use 15 of 16 bytes)

        // this means that v = (uint8_t*)(out + len) - w = (uint8_t*)(out + len) - 1


        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        // pointers
        uint8_t * v = reinterpret_cast<uint8_t *>(out + len) - 1;
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);


        bitgroup_ops<3, BIT_REV_SSSE3, false> op128;

        // since we are working with 120 bits, we are always aligned to offset == 0.
        for (size_t i = 0; i < ((bytes - 1) / 15); i += 15) {
          // enough bytes.  do an iteration
          v -= 15;
          _mm_storeu_si128(reinterpret_cast<__m128i *>(v),
                           op128.reverse<0>(_mm_loadu_si128(reinterpret_cast<__m128i const *>(u))));
          u += 15;
        }


        // take care of remaining bytes.
        if (((bytes + 14) % 15) > 0) {   // bytes - 1 + 15

          // do another iteration with all the remaining.
          if (bytes >= sizeof(__m128i)) {
            // we can just keep 1 completed byte from previous reverse call.
            uint8_t tmp = *(reinterpret_cast<uint8_t *>(out) + sizeof(__m128i) - 1); // get an internal byte from prev reverse.
            _mm_storeu_si128(reinterpret_cast<__m128i *>(out),
                op128.reverse<(2 - ((bytes - 15) * 8) % 3)>(_mm_loadu_si128(reinterpret_cast<__m128i const *>(in + len) - 1)));
            *(reinterpret_cast<uint8_t *>(out) + sizeof(__m128i) - 1) |= tmp;
          } else {
            // too short
            __m128i tmp = op128.reverse<0>(_mm_loadu_si128(reinterpret_cast<__m128i const *>(in)));
            memcpy(out, reinterpret_cast<uint8_t *>(&tmp) + (sizeof(__m128i) - bytes), bytes);
          }

        }  // else the last byte is already taken cared of via loop.
        return (len * sizeof(WORD_TYPE) * 8) % 3;

      }
#else
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint8_t bit_offset = 0 ) {
        //printf("no ssse3. ");
        return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE *>(out),
                                                     ::std::forward<WORD_TYPE const *>(in),
                                                      len,
                                                      bit_offset);

      }
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE, size_t len,
          unsigned int WordsInM128 = sizeof(__m128i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_SSSE3), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
        //printf("no ssse3. ");
        return reverse<BIT_GROUP_SIZE, BIT_REV_SWAR>(::std::forward<WORD_TYPE (&)[len]>(out),
                                                     ::std::forward<WORD_TYPE const (&)[len]>(in));

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
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint8_t bit_offset = 0 ) {

        if ((sizeof(WORD_TYPE) * len) < sizeof(__m256i)) {
          // SSSE3 right now is not performant for power of 2.
          return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>( ::std::forward<WORD_TYPE *>(out),
                                                         ::std::forward<WORD_TYPE const *>(in),
                                                          len,
                                                          bit_offset);
        }


        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= sizeof(__m128i) * 8, "ERROR: currenly reverse does not support 256 BIT_GRUOP_SIZE for SSSE3");



        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
        assert(bit_offset == 0);
//        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) == 0);

        //printf("avx2: ");

        //memset(out, 0, len);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        __m256i * v = reinterpret_cast<__m256i *>(out + len);
        __m256i const * u = reinterpret_cast<__m256i const *>(in);

        //printf("remainder %lu\n", rem);

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2> op256;


        for (; rem >= WordsInM256; rem -= WordsInM256) {
          // enough bytes.  do an iteration
          --v;
          _mm256_storeu_si256(v, op256.reverse(_mm256_loadu_si256(u)));
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
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE,
          typename WORD_TYPE, size_t len,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2) &&
                                          ((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {

        if ((sizeof(WORD_TYPE) * len) < sizeof(__m256i)) {
          // SSSE3 right now is not performant for power of 2.
          return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>( ::std::forward<WORD_TYPE (&)[len]>(out),
                                                         ::std::forward<WORD_TYPE const (&)[len]>(in));
        }


        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= sizeof(__m128i) * 8, "ERROR: currenly reverse does not support 256 BIT_GRUOP_SIZE for SSSE3");



        //if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE requires bit_offset to be 0.");
 //        if (((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
//          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");
        assert(((len * sizeof(WORD_TYPE)) % ((BIT_GROUP_SIZE + 7) / 8)) == 0);

        //memset(out, 0, len);  // needed because we bitwise OR.
        // decide which one to use.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        // pointers
        __m256i * v = reinterpret_cast<__m256i *>(out + len);
        __m256i const * u = reinterpret_cast<__m256i const *>(in);

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2> op256;


        for (size_t i = 0; i < (len / WordsInM256) ;  ++i) {
          // enough bytes.  do an iteration
          --v;
          _mm256_storeu_si256(v, op256.reverse(_mm256_loadu_si256(u)));
          ++u;
        }
        constexpr size_t rem = (len % WordsInM256);
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

        return 0;  // return remainder.
      }
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2) &&
                                          (BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint8_t bit_offset = 0 ) {
        if ((sizeof(WORD_TYPE) * len) < sizeof(__m256i)) {
          // SSSE3 right now is most performant for bit group size of 3
          return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>( ::std::forward<WORD_TYPE *>(out),
                                                         ::std::forward<WORD_TYPE const *>(in),
                                                          len,
                                                          bit_offset);
        }
        size_t bytes = len * sizeof(WORD_TYPE);


        memset(out, 0, bytes);  // needed because we bitwise OR.

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = bytes;
        // pointers
        uint8_t * v = reinterpret_cast<uint8_t *>(out + len);
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);
        uint8_t init_offset = bit_offset;


        bitgroup_ops<3, BIT_REV_AVX2, false> op256;

        for (; rem >= 32; rem -= 32) {
          // enough bytes.  do an iteration
          v -= 32;
          bit_offset = op256.reverse(v, u, 32, bit_offset);
          u += 32;

          if ((rem > 32) && (bit_offset > 0) ) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
        }
        if (rem > 0) {
          // do another iteration with all the remaining.   not that offset is tied to the current u position, so can't use memcpy-avoiding approach.
          if (bytes >= 32) {
            bit_offset = ((bytes - 32) * 8 - init_offset) % BIT_GROUP_SIZE;
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
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE,
        typename WORD_TYPE, size_t len,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2) &&
                                          (BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
        if ((sizeof(WORD_TYPE) * len) < sizeof(__m256i)) {
          // SSSE3 right now is most performant for bit group size of 3
          return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>(::std::forward<WORD_TYPE (&)[len]>(out),
                                                        ::std::forward<WORD_TYPE const (&)[len]>(in));
        }
        constexpr size_t bytes = len * sizeof(WORD_TYPE);


        //memset(out, 0, bytes);  // needed because we bitwise OR.

        // want to unroll a few iterations (y) (small is better for smaller data),
        // so that there is minimal discarded bytes (w) and small overlap bytes (z)
        // (3x + 8w) = ((32-z)y + z) * 8; choose y is 1, z = 2, w = 2.  (i.e, use 30 of 32 bytes)
        // also can choose y = 2, z = 1, w = 0.  logic is more complicated, requires more data
        // choose y = 1.

        // this means that v = (uint8_t*)(out + len) - w = (uint8_t*)(out + len) - 2


        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        // pointers
        uint8_t * v = reinterpret_cast<uint8_t *>(out + len) - 2;
        uint8_t const * u = reinterpret_cast<uint8_t const *>(in);

        bitgroup_ops<3, BIT_REV_AVX2, false> op256;

        for (size_t i = 0; i < ((bytes - 2) / 30); i += 30) {
          // enough bytes.  do an iteration
          v -= 30;
          _mm256_storeu_si256(reinterpret_cast<__m256i *>(v),
                           op256.reverse<0>(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(u))));
          u += 30;
        }

        if (((bytes + 28) % 30) > 0) { // bytes - 2 + 30
          // do another iteration with all the remaining.
          if (bytes >= sizeof(__m256i)) {
            // we can just keep 1 completed byte from previous reverse call.
            uint8_t tmp = *(reinterpret_cast<uint8_t *>(out) + sizeof(__m256i) - 1); // get an internal byte from prev reverse.
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(out),
                op256.reverse<(2 - ((bytes - 31) * 8) % 3)>(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(in + len) - 1)));
            *(reinterpret_cast<uint8_t *>(out) + sizeof(__m256i) - 1) |= tmp;
          } else {
            // too short
            __m256i tmp = op256.reverse<0>(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(in)));
            memcpy(out, reinterpret_cast<uint8_t *>(&tmp) + (sizeof(__m256i) - bytes), bytes);
          }
        }
        return (len * sizeof(WORD_TYPE) * 8) % 3;
      }


#else
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2), uint8_t>::type
      reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len, uint8_t bit_offset = 0 ) {
        //printf("no avx2. ");

        // cascade to SSSE3 and let the choice of implementation be decided there.
        return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>( ::std::forward<WORD_TYPE *>(out),
                                                       ::std::forward<WORD_TYPE const *>(in),
                                                        len,
                                                        bit_offset);
      }
      template <unsigned int BIT_GROUP_SIZE, unsigned char MAX_SIMD_TYPE, typename WORD_TYPE, size_t len,
          unsigned int WordsInM256 = sizeof(__m256i) / sizeof(WORD_TYPE)>
      BITS_INLINE typename std::enable_if<(MAX_SIMD_TYPE == BIT_REV_AVX2), uint8_t>::type
      reverse(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {
        //printf("no ssse3. ");
        return reverse<BIT_GROUP_SIZE, BIT_REV_SSSE3>(::std::forward<WORD_TYPE (&)[len]>(out),
                                                      ::std::forward<WORD_TYPE const (&)[len]>(in));

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
      negate(WORD_TYPE * out, WORD_TYPE const * in, size_t const & len) {

        constexpr unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE);

        // process uint64_t
        size_t len_64 = len - len % WordsInUint64;
        for (size_t i = 0; i < len_64; i += WordsInUint64) {
          *(reinterpret_cast<uint64_t *>(out + i)) = ~(*(reinterpret_cast<uint64_t const *>(in + i)));
        }
        // process remainder.
        for (size_t i = len_64; i < len; ++i) {
          *(out + i) = ~(*(in + i));
        }
      }
      template <typename WORD_TYPE, size_t len>
      BITS_INLINE void
      negate(WORD_TYPE (&out)[len], WORD_TYPE const (&in)[len]) {

        constexpr unsigned int WordsInUint64 = sizeof(uint64_t) / sizeof(WORD_TYPE);

        // process uint64_t
        for (size_t i = 0; i < len / WordsInUint64; ++i) {
          *(reinterpret_cast<uint64_t *>(out) + i) = ~(*(reinterpret_cast<uint64_t const *>(in) + i));
        }
        // process remainder.
        for (size_t i = len - len % WordsInUint64; i < len; ++i) {
          *(out + i) = ~(*(in + i));
        }
      }



    } // namespace bit_ops


  } // namespace utils

} // namespace bliss








#endif /* SRC_UTILS_BITGROUP_OPS_HPP_ */
