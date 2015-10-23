/**
 * @file    bitgroup_ops.hpp
 * @ingroup
 * @author  tpan
 * @brief   collection of bit reverse routines for single word or a byte array (or array of words)
 * @details support reverse by bit groups, i.e. within bit group, the bit order does not change.
 *          Designed as a set of functors because this allows us to do partial template specialization (BIT GROUP SIZE and WORD_TYPE, whereas functions need to be fully specialized.
 *
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */
#ifndef SRC_UTILS_BITGROUP_OPS_HPP_
#define SRC_UTILS_BITGROUP_OPS_HPP_

#include <type_traits>   // is_integral, etc
#include <limits>        // numeric_limits
#include <string>
#include <cstring>       // memcpy and memset
#include <sstream>  // stringstream
#include <cstdint>       // uint8_t

#include "utils/logging.h"

#include <x86intrin.h>   // all intrinsics.  will be enabled based on compiler flag such as __SSSE3__ internally.


// done:  3bit reverse - done for SWAR, SSSE3, AVX2
// done:  see effect of not shifting when working with byte arrays.  okay not to shift, if bit_offset is passed to reverse
// done:  vector reverse.
// TODO:  perf compare to old impl.
// TODO:  add other operations such as masking, shifting, loading data from mem and storing data to mem.  refactor array version of reverse to use load/op/store


namespace bliss {

  namespace utils {

    namespace bit_ops {

      static constexpr unsigned char BIT_REV_SEQ = 0;
      static constexpr unsigned char BIT_REV_SWAR = 1;   // SIMD Within A Register

#if defined(__SSSE3__)
      static constexpr unsigned char BIT_REV_SSSE3 = 2;
#endif

#if defined(__AVX2__)
      static constexpr unsigned char BIT_REV_AVX2 = 4;
#endif


      /// shift right by number of bits less than 8.
      template <typename WORD_TYPE>
      inline WORD_TYPE srli(WORD_TYPE const val, uint8_t shift) {
        return val >> shift;
      }
      /// shift left by number of bits less than 8.
      template <typename WORD_TYPE>
      inline WORD_TYPE slli(WORD_TYPE const val, uint8_t shift) {
        return val << shift;
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
          static_assert(BIT_GROUP_SIZE <= (sizeof(size_t) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in size_t");
          static_assert(((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0) || (BIT_GROUP_SIZE < 8), "ERROR: BIT_GROUP_SIZE is has to be powers of 2 or less than 8");

          /* Linear (inefficient) reverse: */

          static constexpr size_t group_mask = ~(std::numeric_limits<size_t>::max() << BIT_GROUP_SIZE);

          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_SIMD_TYPE;

          /**
           * @brief reverse function to reverse bits for data types are are not power of 2 number of bytes
           * @details OR the current len size result with any existing output, so important to make sure the underlying out vector is appropriately initialized.
           *          len has to be a multiple of BIT_GROUP_SIZE if BIT_GROUP_SIZE is greater than 8 (i.e. multiple bytes).
           *          in other words, when a converted word is divided into BIT_GROUP_SIZE, there should not be padding zero's in the groups
           *
           *
           * @param out       pointer into a byte array.  the array should be initialized to 0 before any call to this function.
           * @return number of bits in the last byte that are not swapped (i.e. last 2 bits for 3 bits op.)
           */
          size_t reverse(uint8_t * out, uint8_t const * in, size_t const len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be less than 8.  for larger, shift input instead.");

            if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
              throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");

//            std::cout << "in: ";
//            for (size_t i = 0; i < len; ++i) {
//              std::cout << std::hex << static_cast<size_t>(in[len - 1 - i]) << " ";
//
//            }
//            std::cout << std::endl;

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
                  uint64_t v7 = 0;
                  memcpy(reinterpret_cast<uint8_t*>(&v7), in, len);
                  ///v7 >>= bit_offset;
                  v7 = this->reverse(v7, bit_offset); // >> bit_offset;
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
                  uint32_t v3 = 0;
                  memcpy(reinterpret_cast<uint8_t*>(&v3), in, 3);
                  // v3 >>= bit_offset;
                  v3 = this->reverse(v3, bit_offset); // >> bit_offset;
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
                ERRORF("ERROR:  length of array should be > 0 and <= 8.  For longer, break it up or use the vector version of bitgroup_ops.");
                break;
            }
//            std::cout << "out: ";
//            for (size_t i = 0; i < len; ++i) {
//              std::cout << std::hex << static_cast<size_t>(out[len - 1 - i]) << " ";
//
//            }
//            std::cout << std::endl;

            return (len * 8 - bit_offset) % BIT_GROUP_SIZE;
          }

          /**
           * @brief  reverse the lower portion of WORD, starting from 0, up to the largest multiples of BIT_GROUP_SIZE that still fits.
           *
           * @param v    input word.  LSB aligned to start of a bit group.
           * @return     bit reversed word.  MSB aligned to end of the originally lowest bit group.
           */
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS >= (sizeof(WORD_TYPE) * 8)), int>::type = 1>
          WORD_TYPE reverse(WORD_TYPE const & v, uint8_t bit_offset = 0) {
//            static_assert(BIT_GROUP_SIZE == (sizeof(WORD_TYPE) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");
            return v;
          }
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS < (sizeof(WORD_TYPE) * 8)), int>::type = 1>
          WORD_TYPE reverse(WORD_TYPE const & v, uint8_t bit_offset = 0) {

            // copy the source data
            WORD_TYPE w = 0;
            WORD_TYPE u = v;

            unsigned int rem = (sizeof(WORD_TYPE) * 8 - bit_offset) % BIT_GROUP_SIZE;
            //WORD_TYPE rem_mask = static_cast<WORD_TYPE>(~(std::numeric_limits<WORD_TYPE>::max() >> rem));
            //std::cout << "begin in: " << std::hex << static_cast<size_t>(u) << " out: " << static_cast<size_t>(w) << std::endl;

            u >>= bit_offset;

            // now do the middle part
            for (size_t i = bit_offset; i < (sizeof(WORD_TYPE) * 8 - rem); i += BIT_GROUP_SIZE) {
              w <<= BIT_GROUP_SIZE;
              w |= (u & group_mask);
              u >>= BIT_GROUP_SIZE;

              //std::cout << " in: " << std::hex << static_cast<size_t>(u) << " out: " << static_cast<size_t>(w) << std::endl;
            }

            //std::cout << "rem mask: " << rem_mask << std::endl;

            // don't put back the remainder.  just shift.. we'll need to OR bytes together anyways.
            w <<= rem;
//            w |= (u & rem_mask);
            //std::cout << "end in: " << std::hex << static_cast<size_t>(u) << " out: " << static_cast<size_t>(w) << std::endl;

            // finally merge back to v.
            return w;
          }


      };

      /// partial template specialization for SWAR based bit reverse.  uses BSWAP when appropriate, else use bit shift  this is defined only for bit_group_sizes that are powers of 2 up to sizeof WORD_TYPE)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SWAR, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE <= (sizeof(size_t) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in size_t");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");


          // masks.  if size_t indicates that system is 32 bit, then this should truncate.
          static constexpr size_t mask32 = 0x00000000FFFFFFFF;
          static constexpr size_t mask16 = 0x0000FFFF0000FFFF;
          static constexpr size_t  mask8 = 0x00FF00FF00FF00FF;
          static constexpr size_t  mask4 = 0x0F0F0F0F0F0F0F0F;
          static constexpr size_t  mask2 = 0x3333333333333333;
          static constexpr size_t  mask1 = 0x5555555555555555;

          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_SWAR;

          /// reverse function to reverse bits for data types are are not power of 2 number of bytes
          size_t reverse(uint8_t * out, uint8_t const * in, size_t len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            if (bit_offset > 0) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be 0 for SWAR with power of 2 BIT_GROUP_SIZE.");


            if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
              throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");

//            std::cout << "in: ";
//            for (size_t i = 0; i < len; ++i) {
//              std::cout << std::hex << static_cast<size_t>(in[len - 1 - i]) << " ";
//
//            }
//            std::cout << std::endl;
            switch (len) {
              case 8:
                *(reinterpret_cast<uint64_t*>(out)) = this->reverse(*(reinterpret_cast<const uint64_t*>(in)));
                break;
              case 7:
              case 6:
              case 5:
                {
                  uint64_t v64 = 0;
                  memcpy(reinterpret_cast<uint8_t*>(&v64), in, len);
                  v64 = this->reverse(v64);
                  memcpy(out, reinterpret_cast<uint8_t*>(&v64) + (8-len), len);
                }
                break;
              case 4:
                *(reinterpret_cast<uint32_t*>(out)) = this->reverse(*(reinterpret_cast<const uint32_t*>(in)));
                break;
              case 3:
                {
                  uint32_t v32 = 0;
                  memcpy(reinterpret_cast<uint8_t*>(&v32), in, len);
                  v32 = this->reverse(v32);
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
                ERRORF("ERROR:  length of array should be > 0 and <= 8.  For longer, break it up or use the vector version of bitgroup_ops.");
                break;

            }

//            std::cout << "out: ";
//            for (size_t i = 0; i < len; ++i) {
//              std::cout << std::hex << static_cast<size_t>(out[len - 1 - i]) << " ";
//
//            }
//            std::cout << std::endl;

            return 0;
          }

          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS >= (sizeof(WORD_TYPE) * 8)), int>::type = 1>
          WORD_TYPE reverse(WORD_TYPE const &u) {
//            static_assert(BIT_GROUP_SIZE == (sizeof(WORD_TYPE) * 8), "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }
          template <typename WORD_TYPE, unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS < (sizeof(WORD_TYPE) * 8)), int>::type = 1>
          WORD_TYPE reverse(WORD_TYPE const &u) {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            WORD_TYPE v = u;

            // essentially a Duff's Device here, cascading to smaller swap sizes.
            // swap the bit_groups in logarithmic steps.  only execute steps that have swap size smaller than the WORD_TYPE, and only if bit group is <= the current swap size
            if (BIT_GROUP_SIZE <= 8) {
              switch (sizeof(WORD_TYPE)) {
                case 8:
                  DEBUGF("BSWAP reverse for 64 bit word");
                  v = _bswap64(v);
                  break;    // swap all 8 bytes
                case 4:
                  DEBUGF("BSWAP reverse for 32 bit word");
                  v = _bswap(v);
                  break;  // swap all 4 bytes
                case 2:
                  DEBUGF("BSWAP reverse for 16 bit word");
                  v = ((v >>  8) &  mask8) | ((v &  mask8) <<  8);
                  break;  // swap both 2 byte
              }
            } else {  // BIT_GROUP_SIZE > 8:  16 or 32, depending on word size..  limited by word_type size
              switch (BIT_GROUP_SIZE) {
                // should be shifting by 32 but causes compiler warning if v is same size of BIT_GROUP_SIZE
                case 16:
                  switch (sizeof(WORD_TYPE)) { // instead of 32 and 16, do BIT_GROUP_SIZE 2x
                    case 8:
                      DEBUGF("SWAR reverse for in 32 group");
                      v = (((v >> BIT_GROUP_SIZE) >> BIT_GROUP_SIZE) & mask32) | (((v & mask32) << BIT_GROUP_SIZE) << BIT_GROUP_SIZE );         // swap 4 bytes
                    case 4:
                      DEBUGF("SWAR reverse in 16 bit groups");
                      v = ((v >> BIT_GROUP_SIZE) & mask16) | ((v & mask16) << BIT_GROUP_SIZE);
                      break;  // swap 2 bytes
                  };
                  break;
                    case 32:
                      DEBUGF("SWAR reverse in 32 bit groups");
                      v = ((v >> BIT_GROUP_SIZE) & mask32) | ((v & mask32) << BIT_GROUP_SIZE);
                      break;  // swap 4 bytes.  wordtype has to be size_t
              }
            }

            // now swap if requested group size is 4, 2, or 1.  only execute if bit group is <= to current swap size.
            switch (BIT_GROUP_SIZE) {
              case 1:
                DEBUGF("SWAR reverse in 1 bit groups");
                v = ((v >>  1) &  mask1) | ((v &  mask1) <<  1);  // swap 1 bits
              case 2:
                DEBUGF("SWAR reverse in 2 bit groups");
                v = ((v >>  2) &  mask2) | ((v &  mask2) <<  2);  // swap 2 bits
              case 4:
                DEBUGF("SWAR reverse in 4 bit groups");
                v = ((v >>  4) &  mask4) | ((v &  mask4) <<  4);  // swap nibbles
            }

            return v;
          }

      };


      /// full template specialization for SWAR based bit reverse for 3 bits.  uses BSWAP when appropriate, else use bit shift.  this performs a SWAR style bit shuffling after a full bitwise reverse.)
      template <bool POW2>
      struct bitgroup_ops<3, BIT_REV_SWAR, POW2> {

          // masks.  if size_t indicates that system is 32 bit, then this should truncate.
          static constexpr size_t  mask8 = 0x00FF00FF00FF00FF;
          static constexpr size_t  mask4 = 0x0F0F0F0F0F0F0F0F;
          static constexpr size_t  mask2 = 0x3333333333333333;
          static constexpr size_t  mask1 = 0x5555555555555555;

          static constexpr size_t mask3lo = 0x9249249249249249;
          static constexpr size_t mask3mid = 0x2492492492492492;
          static constexpr size_t mask3hi = 0x4924924924924924;

          static constexpr unsigned int bitsPerGroup = 3;
          static constexpr unsigned char simd_type = BIT_REV_SWAR;

          /// reverse function to reverse bits for data types are are not power of 2 number of bytes
          size_t reverse(uint8_t * out, uint8_t const * in, size_t len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be less than 8.  for larger, shift input instead.");

//            std::cout << "in: ";
//            for (size_t i = 0; i < len; ++i) {
//              std::cout << std::hex << static_cast<size_t>(in[len - 1 - i]) << " ";
//
//            }
//            std::cout << std::endl;

            // SHIFTS are NOT NECESSARY.



            switch (len) {
              case 8:
                {
                  uint64_t v8 = *(reinterpret_cast<const uint64_t*>(in));
                  //v8 >>= bit_offset;
                  *(reinterpret_cast<uint64_t*>(out)) |= this->reverse(v8, bit_offset); // >> bit_offset;
                }
                break;
              case 7:
              case 6:
              case 5:
                {
                  uint64_t v7 = 0;
                  memcpy(reinterpret_cast<uint8_t*>(&v7), in, len);
                  //v7 >>= bit_offset;
                  v7 = this->reverse(v7, bit_offset); // >> bit_offset;
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
                  uint32_t v3 = 0;
                  memcpy(reinterpret_cast<uint8_t*>(&v3), in, 3);
                  //v3 >>= bit_offset;
                  v3 = this->reverse(v3, bit_offset); // >> bit_offset;
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
                  //v1 >>= bit_offset;
                  *out |= this->reverse(v1, bit_offset); // >> bit_offset;
                }
                break;
              default:
                ERRORF("ERROR:  length of array should be > 0 and <= 8.  For longer, break it up or use the vector version of bitgroup_ops.");
                break;
            }

//            std::cout << "out: ";
//            for (size_t i = 0; i < len; ++i) {
//              std::cout << std::hex << static_cast<size_t>(out[len - 1 - i]) << " ";
//
//            }
//            std::cout << std::endl;

            return (len * 8 - bit_offset) % 3;
          }

          template <typename WORD_TYPE>
          WORD_TYPE reverse(WORD_TYPE const &u, uint8_t bit_offset = 0) {
            static_assert((::std::is_integral<WORD_TYPE>::value) && (!::std::is_signed<WORD_TYPE>::value), "ERROR: WORD_TYPE has to be unsigned integral type.");
            static_assert(sizeof(WORD_TYPE) <= 8, "ERROR: WORD_TYPE should be primitive and smaller than 8 bytes for SWAR");

            WORD_TYPE v = u;

            // essentially a Duff's Device here, cascading to smaller swap sizes.
            // swap the bit_groups in logarithmic steps.  only execute steps that have swap size smaller than the WORD_TYPE, and only if bit group is <= the current swap size
            switch (sizeof(WORD_TYPE)) {
              case 8:
                DEBUGF("BSWAP reverse for 64 bit word");
                v = _bswap64(v);
                break;    // swap all 8 bytes
              case 4:
                DEBUGF("BSWAP reverse for 32 bit word");
                v = _bswap(v);
                break;  // swap all 4 bytes
              case 2:
                DEBUGF("BSWAP reverse for 16 bit word");
                v = ((v >>  8) &  mask8) | ((v &  mask8) <<  8);
                break;  // swap both 2 byte
            } // else single byte, no byte swap needed.

            // now swap the bits as if requested group size is 1
            DEBUGF("SWAR reverse in 4 bit groups");
            v = ((v >>  4) &  mask4) | ((v &  mask4) <<  4);  // swap nibbles
            DEBUGF("SWAR reverse in 2 bit groups");
            v = ((v >>  2) &  mask2) | ((v &  mask2) <<  2);  // swap 2 bits
            DEBUGF("SWAR reverse in 1 bit groups");
            v = ((v >>  1) &  mask1) | ((v &  mask1) <<  1);  // swap 1 bits

            // and then swap the 1st and 3rd elements in the group.  this should suffice.
            // note that remainders are at low bits now, so we need 3 patterns for shifting.

            // 2 shifts, 2 ors, 3 ands
            if (((sizeof(WORD_TYPE) * 8 - bit_offset) % 3) == 2) {
              // rem == 2:                    mid, hi, lo
              v = ((v &  mask3mid) >>  2) | (v & mask3lo ) | ((v &  mask3hi ) << 2);
            } else if (((sizeof(WORD_TYPE) * 8 - bit_offset) % 3) == 1) {
              // rem == 1:                    hi, lo, mid
              v = ((v &  mask3lo ) >>  2) | (v & mask3hi ) | ((v &  mask3mid) << 2);
            } else {
              // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
              v = ((v &  mask3hi ) >>  2) | (v & mask3mid) | ((v &  mask3lo ) << 2);
            }


///            std::cout << std::hex << "bit offset = " << static_cast<size_t>(bit_offset) << " in "  << static_cast<size_t>(u) << " out " << static_cast<size_t>(v) << std::endl;

            return v;
          }

      };

#if defined(__SSSE3__)



      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 128bit)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE <= 128, "ERROR: BIT_GROUP_SIZE is greater than number of bits in __m128i");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");



          /// index for reversing the bytes in groups of varying size (8: 1 byte; 16: 2 bytes, etc.)
          static constexpr uint8_t  rev_idx8[16] alignas(16) = {0x0F,0x0E,0x0D,0x0C,0x0B,0x0A,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,0x00};
          static constexpr uint8_t rev_idx16[16] alignas(16) = {0x0E,0x0F,0x0C,0x0D,0x0A,0x0B,0x08,0x09,0x06,0x07,0x04,0x05,0x02,0x03,0x00,0x01};
          static constexpr uint8_t rev_idx32[16] alignas(16) = {0x0C,0x0D,0x0E,0x0F,0x08,0x09,0x0A,0x0B,0x04,0x05,0x06,0x07,0x00,0x01,0x02,0x03};
          static constexpr uint8_t rev_idx64[16] alignas(16) = {0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07};

          /// mask for lower 4 bits
          static constexpr uint8_t mask4[16] alignas(16) = {0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F};

          /// lookup table for reversing bits in a byte in bit groups of 1
          static constexpr uint8_t lut1_lo[16] alignas(16) = {0x00,0x08,0x04,0x0c,0x02,0x0a,0x06,0x0e,0x01,0x09,0x05,0x0d,0x03,0x0b,0x07,0x0f};
          static constexpr uint8_t lut1_hi[16] alignas(16) = {0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0};

          /// lookup table for reversing bits in a byte in bit groups of 2
          static constexpr uint8_t lut2_lo[16] alignas(16) = {0x00,0x04,0x08,0x0c,0x01,0x05,0x09,0x0d,0x02,0x06,0x0a,0x0e,0x03,0x07,0x0b,0x0f};
          static constexpr uint8_t lut2_hi[16] alignas(16) = {0x00,0x40,0x80,0xc0,0x10,0x50,0x90,0xd0,0x20,0x60,0xa0,0xe0,0x30,0x70,0xb0,0xf0};

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_SSSE3;

          static ::std::string toString(__m128i const & v) {
            uint8_t tmp[16] alignas(16);
            _mm_store_si128((__m128i*)tmp, v);
            ::std::stringstream ss;
            for (int i = 15; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }

          /// reverse function to reverse bits for data types are are not __mm128i  (has to be bigger than at least half)
          size_t reverse(uint8_t * out, uint8_t const * in, size_t len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
              throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");

            if (len > 16) throw std::invalid_argument("ERROR:  length of array should be <= 16.  For longer, use a different bitgroup_ops (AVX2) or break it up.");
            //assert(len > 8); //  "ERROR:  length of array should be > 8.  For shorter, use a different bitgroup_ops type.");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            if (bit_offset > 0) throw ::std::invalid_argument("ERROR reversing byte array:  bit offset should be 0 for SSSE3 bit reverse and power of 2 BIT_GROUP_SIZE.");


            if (len == 16) {
            	_mm_storeu_si128((__m128i*)out, this->reverse( _mm_loadu_si128((__m128i*)in) ));
            } else {
				// convert to __m128i, then call the __m128i version.  note that there aren't types that have 128 bits other than simd types
				// allocate temp space
				uint8_t tmp[16] alignas(16);

				// copy data to end of array (reverse will move it to front.
				memcpy(tmp, in, len);

				// copy in, do reverse, copy back out
				_mm_store_si128((__m128i*)tmp, this->reverse( _mm_load_si128((__m128i*)tmp) ));                    // SSE2

				// maskmoveu is very slow.  why?
				//_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2
				memcpy(out, tmp + 16 - len, len);
            }
            return 0;
          }


          template <unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS >= 128), int>::type = 1>
          __m128i reverse(__m128i const & u) {
            static_assert(BIT_GROUP_SIZE == 128, "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }

          // groupsize =1 version:  4 loads, 3 shuffles, 3 logical ops, 1 shift op
          template <unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS < 128), int>::type = 1>
          __m128i reverse(__m128i const & u) {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            //        print(v, "v");

            __m128i v = u;

            //== first reverse the bit groups via 1 shuffle operation
            switch (BIT_GROUP_SIZE) {
              case 1:
              case 2:
              case 4:
              case 8:
                // if bit group is <= 8, then need to shuffle bytes, and then shuffle bits.
                DEBUGF("SSSE3 reverse in 8 bit groups");
                v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*) rev_idx8));
                break;                               // SSSE3
              case 16:   // if bit group is 16, then shuffle using epi8
                DEBUGF("SSSE3 reverse in 16 bit groups");
                v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*)rev_idx16));
                break;                                // SSSE3
              case 32: //v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*)rev_idx32));                               // SSSE3
                DEBUGF("SSSE3 reverse in 32 bit groups");
                v = _mm_shuffle_epi32(v, 0x1B);
                break; // original 11 10 01 00 (high to low) => 00 01 10 11 (high to low) == 0x1B                // SSE2
              case 64:  //v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*)rev_idx64));                               // SSSE3
                DEBUGF("SSSE3 reverse in 64 bit groups");
                v = _mm_shuffle_epi32(v, 0x4E);
                break;  // original 11 10 01 00 (high to low) => 01 00 11 10 (high to low) == 0x4E				   // SSE2
            }
            //        print(r, "r");

            //== now see if we need to shuffle bits.
            if (BIT_GROUP_SIZE < 8) {

              //== next get the lut indices.
              // load constants:
              __m128i _mask_lo = _mm_load_si128((__m128i*)mask4);                                         // SSE2

              //        print(mask_lo, "mask_lo");

              // lower 4 bits
              __m128i lo = _mm_and_si128(_mask_lo, v);                                                      // SSE2
              // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
              // shift not using si128 - that instruction shifts bytes not bits.
              __m128i hi = _mm_srli_epi16(_mm_andnot_si128(_mask_lo, v), 4);                                // SSE2

              //        print(lo, "lo");
              //        print(hi, "hi");

              switch (BIT_GROUP_SIZE) {
                case 1:
                  DEBUGF("SSSE3 reverse in 1 bit groups");

                  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                  lo = _mm_shuffle_epi8(_mm_load_si128((__m128i*)lut1_hi), lo);                        // SSSE3
                  hi = _mm_shuffle_epi8(_mm_load_si128((__m128i*)lut1_lo), hi);                        // SSSE3
                  break;
                case 2:
                  DEBUGF("SSSE3 reverse in 2 bit groups");
                  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                  lo = _mm_shuffle_epi8(_mm_load_si128((__m128i*)lut2_hi), lo);                        // SSSE3
                  hi = _mm_shuffle_epi8(_mm_load_si128((__m128i*)lut2_lo), hi);                        // SSSE3
                  break;
                case 4:
                  DEBUGF("SSSE3 reverse in 4 bit groups");
                  lo = _mm_slli_epi16(lo, 4);  // shift lower 4 to upper
                  break;
              }

              // recombine
              v =  _mm_or_si128(lo, hi);                                                                  // SSE2

              //      print(res, "result");
            }

            return v;
          }


      };

      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::rev_idx8[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::rev_idx16[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::rev_idx32[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::rev_idx64[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::mask4[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::lut1_lo[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::lut1_hi[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::lut2_lo[16];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3, POW2>::lut2_hi[16];


      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 128bit)
      template <bool POW2>
      struct bitgroup_ops<3, BIT_REV_SSSE3, POW2> {

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = 3;
          static constexpr unsigned char simd_type = BIT_REV_SSSE3;

          static constexpr size_t mask3lo[2] alignas(16) = { 0x9249249249249249, 0x4924924924924924 };
          static constexpr size_t mask3mid[2] alignas(16) = { 0x2492492492492492, 0x9249249249249249 };
          static constexpr size_t mask3hi[2] alignas(16) = { 0x4924924924924924, 0x2492492492492492 };

          bitgroup_ops<1, BIT_REV_SSSE3, true> bit_rev_1;

          static ::std::string toString(__m128i const & v) {
            uint8_t tmp[16] alignas(16);
            _mm_store_si128((__m128i*)tmp, v);
            ::std::stringstream ss;
            for (int i = 15; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }

          /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
          __m128i srli(__m128i val, uint8_t shift) {
            // get the 9th byte into the 8th position
            // shift by 1 byte.  avoids a mask load during shuffle, and is friendly with epi16, epi32, and epi64 versions of shifts.
            __m128i tmp = _mm_srli_si128(val, 1);

            // then left shift to get the bits in the right place
            // shift the input value via epi64 by the number of shifts
            // then or together.
            return _mm_or_si128(_mm_slli_epi64(tmp, 8 - shift), _mm_srli_epi64(val, shift));
          }

          /// shift left by number of bits less than 8.
          __m128i slli(__m128i val, uint8_t shift) {
            // get the 9th byte into the 8th position
            // shift by 1 byte.  avoids a mask load during shuffle, and is friendly with epi16, epi32, and epi64 versions of shifts.
            __m128i tmp = _mm_slli_si128(val, 1);

            // then left shift to get the bits in the right place
            // shift the input value via epi64 by the number of shifts
            // then or together.
            return _mm_or_si128(_mm_srli_epi64(tmp, 8 - shift), _mm_slli_epi64(val, shift));
          }


          /// reverse function to reverse bits for data types are are not __mm128i  (has to be bigger than at least half)
          size_t reverse(uint8_t * out, uint8_t const * in, size_t len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            if (len > 16) throw std::invalid_argument("ERROR:  length of array should be <= 16.  For longer, use a different bitgroup_ops (AVX2) or break it up.");
            //assert(len > 8); //  "ERROR:  length of array should be > 8.  For shorter, use a different bitgroup_ops type.");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array SSSE3:  bit offset should be less than 8.  for larger, shift input instead.");

			// convert to __m128i, then call the __m128i version.  note that there aren't types that have 128 bits other than simd types

            // allocate temp space
            uint8_t tmp[16] alignas(16);

            __m128i v;
            __m128i w;

            // load the data
            if (len == 16) {
            	v = _mm_loadu_si128((__m128i*)in);
            } else {
            	memcpy(tmp, in, len);

            	v = _mm_load_si128((__m128i*)tmp);
            }

            // do the conversion
            if (bit_offset == 0) {
            	w = this->reverse( v );
            } else {
            	// ARE THE SHIFTS NECESSARY?  they end up zeroing out the offset, which could be accomplished via a mask.
            	// note that the remainder portion is not zeroed, specifically, if rem == 2, then the center bit is not zeroed, but the side bits are zeroed.
              // center bit OR with itself does not change.  side bits from next iter ORed with zeroed side bits (this iter) returns the right results.
            	// SO ARE THE ZEROING ALSO NECESSARY?

              // SHIFTS are NOT NECESSARY.



            	// shift right before
//            	v = srli(v, bit_offset);

            	// convert
            	w = this->reverse( v , bit_offset);

            	// shift right after
//            	w = srli(w, bit_offset);
            }

            // copy back to out.  ALWAYS NEED TO OR, so can never directly store.
            _mm_store_si128((__m128i*)tmp, w);                    // SSE2
            if (len > 2) memcpy(out+1, tmp + (17 - len), len - 2);
            if (len > 1) out[len - 1] |= tmp[15];
            out[0] |= tmp[16 - len];

            // maskmoveu is very slow.  why?
			//_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2

            // return remainder.
            return (len * 8 - bit_offset) % 3;
          }


          __m128i reverse(__m128i const & u, uint8_t bit_offset = 0 ) {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            //        print(v, "v");

            // first reverse bits in groups of 1.
            __m128i v = bit_rev_1.reverse(u);  // 4 loads, 3 shuffles, 3 logical ops, 1 shift op

            // then get the individual bits.  3 loads, 3 ORs.
            __m128i lo = _mm_and_si128(v, _mm_load_si128((__m128i*)mask3lo));
            __m128i mid = _mm_and_si128(v, _mm_load_si128((__m128i*)mask3mid));
            __m128i hi = _mm_and_si128(v, _mm_load_si128((__m128i*)mask3hi));

            // next shift based on the rem bits.  2 cross boundary shifts = 2 byteshift, 4 shifts, 2 ORs. + 2 ORs.
            if (((128 - bit_offset) % 3) == 2) {
              // rem == 2:                    mid, hi, lo
              v = _mm_or_si128(srli(mid, 2), _mm_or_si128(lo, slli( hi, 2) ) );
            } else if (((128 - bit_offset) % 3) == 1) {
              // rem == 1:                    hi, lo, mid
                v = _mm_or_si128(srli(lo, 2), _mm_or_si128(hi, slli( mid, 2) ) );
            } else {
              // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
                v = _mm_or_si128(srli(hi, 2), _mm_or_si128(mid, slli( lo, 2) ) );
            }

            // alternatively,  for cross-boundary shift, do byte shift, and reverse shift (8-bit_offset), then or. with shifted (offset)->  each CBS is finished in 4 ops w/o load.


            return v;
          }


      };

      // need the DUMMY template parameter to correctly instantiate here.
      template <bool POW2>
      constexpr size_t bitgroup_ops<3, BIT_REV_SSSE3, POW2>::mask3lo[2];
      template <bool POW2>
      constexpr size_t bitgroup_ops<3, BIT_REV_SSSE3, POW2>::mask3mid[2];
      template <bool POW2>
      constexpr size_t bitgroup_ops<3, BIT_REV_SSSE3, POW2>::mask3hi[2];



#endif

#if defined(__AVX2__)




      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 256bit)
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      struct bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2> {
          static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE is 0");
          static_assert(BIT_GROUP_SIZE <= 256, "ERROR: BIT_GROUP_SIZE is greater than number of bits in __m256i");
          static_assert((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0, "ERROR: BIT_GROUP_SIZE is has to be powers of 2");


          /// index for reversing the bytes in groups of varying size (8: 1 byte; 16: 2 bytes, etc.)
          static constexpr uint8_t  rev_idx_lane8[32] alignas(32) = {0x0F,0x0E,0x0D,0x0C,0x0B,0x0A,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,0x00,
                                                                     0x0F,0x0E,0x0D,0x0C,0x0B,0x0A,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,0x00};
          static constexpr uint8_t rev_idx_lane16[32] alignas(32) = {0x0E,0x0F,0x0C,0x0D,0x0A,0x0B,0x08,0x09,0x06,0x07,0x04,0x05,0x02,0x03,0x00,0x01,
                                                                     0x0E,0x0F,0x0C,0x0D,0x0A,0x0B,0x08,0x09,0x06,0x07,0x04,0x05,0x02,0x03,0x00,0x01};
          static constexpr uint8_t rev_idx32[32] alignas(32) = {0x07,0x00,0x00,0x00,0x06,0x00,0x00,0x00,0x05,0x00,0x00,0x00,0x04,0x00,0x00,0x00,
                                                                0x03,0x00,0x00,0x00,0x02,0x00,0x00,0x00,0x01,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
          //          static constexpr uint8_t rev_idx64[32] alignas(32) = {0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
          //        		  0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07};


          /// mask for lower 4 bits
          static constexpr uint8_t mask4[32] alignas(32) = {0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,
                                                            0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F};

          /// lookup table for reversing bits in a byte in bit groups of 1
          static constexpr uint8_t lut1_lo[32] alignas(32) = {0x00,0x08,0x04,0x0c,0x02,0x0a,0x06,0x0e,0x01,0x09,0x05,0x0d,0x03,0x0b,0x07,0x0f,
                                                              0x00,0x08,0x04,0x0c,0x02,0x0a,0x06,0x0e,0x01,0x09,0x05,0x0d,0x03,0x0b,0x07,0x0f};
          static constexpr uint8_t lut1_hi[32] alignas(32) = {0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0,
                                                              0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0};

          /// lookup table for reversing bits in a byte in bit groups of 2
          static constexpr uint8_t lut2_lo[32] alignas(32) = {0x00,0x04,0x08,0x0c,0x01,0x05,0x09,0x0d,0x02,0x06,0x0a,0x0e,0x03,0x07,0x0b,0x0f,
                                                              0x00,0x04,0x08,0x0c,0x01,0x05,0x09,0x0d,0x02,0x06,0x0a,0x0e,0x03,0x07,0x0b,0x0f};
          static constexpr uint8_t lut2_hi[32] alignas(32) = {0x00,0x40,0x80,0xc0,0x10,0x50,0x90,0xd0,0x20,0x60,0xa0,0xe0,0x30,0x70,0xb0,0xf0,
                                                              0x00,0x40,0x80,0xc0,0x10,0x50,0x90,0xd0,0x20,0x60,0xa0,0xe0,0x30,0x70,0xb0,0xf0};

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = BIT_GROUP_SIZE;
          static constexpr unsigned char simd_type = BIT_REV_AVX2;

          static ::std::string toString(__m256i const & v) {
            uint8_t tmp[32] alignas(32);
            _mm256_store_si256((__m256i*)tmp, v);																// AVX
            ::std::stringstream ss;
            for (int i = 31; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }


          /// reverse function to reverse bits for data types are are not __mm256i   (has to be bigger than at least half)
          size_t reverse(uint8_t * out, uint8_t const * in, size_t len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
              throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");


            if (len > 32) throw std::invalid_argument("ERROR:  length of array should be <= 32.  For longer, break it up.");

            // convert to __m128i, then call the __m128i version.  note that there aren't types that have 128 bits other than simd types
            // allocate temp space
            uint8_t tmp[32] alignas(32);

            // copy data to end of array (reverse will move it to front.
            memcpy(tmp + 32 - len, in, len);

            // copy in, do reverse, copy back out
            // NOTE: aligned load does not seem to work correctly - SEGV.  use unaligned load.
            _mm256_storeu_si256((__m256i*)tmp, this->reverse( _mm256_loadu_si256((__m256i*)tmp) ));                    // AVX

            // maskmoveu is very slow.  why?
            //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2
            memcpy(out, tmp, len);

            return 0;
          }

          template <unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS >= 256), int>::type = 1>
          __m256i reverse(__m256i const & u) {
            static_assert(BIT_GROUP_SIZE == 256, "ERROR: BIT_GROUP_SIZE is greater than number of bits in WORD_TYPE");

            return u;
          }
          template <unsigned int BITS = BIT_GROUP_SIZE, typename std::enable_if<(BITS < 256), int>::type = 1>
          __m256i reverse(__m256i const & u) {


            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            //        print(v, "v");
            __m256i v = u;

            //== first reverse the bit groups via 1 shuffle operation
            switch (BIT_GROUP_SIZE) {
              case 1:
              case 2:
              case 4:
              case 8:
                // if bit group is <= 8, then need to shuffle bytes, and then shuffle bits.
                DEBUGF("AVX2 reverse in 8 bit groups");
                v = _mm256_shuffle_epi8(v, _mm256_load_si256((__m256i*) rev_idx_lane8));     // reverse 8 in each of 128 bit lane          // AVX2
                v = _mm256_permute2x128_si256(v, v, 0x03);							    // then swap the lanes                           // AVX2.  latency = 3
                break;
              case 16:
                DEBUGF("AVX2 reverse in 16 bit groups");
                v = _mm256_shuffle_epi8(v, _mm256_load_si256((__m256i*)rev_idx_lane16));     // reverse 16 in each of 128 bit lane         // AVX2
                v = _mm256_permute2x128_si256(v, v, 0x03);							    // then swap the lanes							 // AVX2  latency = 3
                break;
              case 32:
                //v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*)rev_idx32));                               // SSSE3
                DEBUGF("AVX2 reverse in 32 bit groups");
                v = _mm256_permutevar8x32_epi32(v, _mm256_load_si256((__m256i*)rev_idx32));                         // AVX2 latency = 3
                break; // original 11 10 01 00 (high to low) => 00 01 10 11 (high to low) == 0x1B          // AVX2
              case 64:
                //v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*)rev_idx64));                               // SSSE3
                DEBUGF("AVX2 reverse in 64 bit groups");
                v = _mm256_permute4x64_epi64(v, 0x1B);                                                              // AVX2 latency = 3
                break; // original 11 10 01 00 (high to low) => 00 01 10 11 (high to low) == 0x1B         // AVX2
              case 128:
                //v = _mm_shuffle_epi8(v, _mm_load_si128((__m128i*)rev_idx64));                               // SSSE3
                DEBUGF("AVX2 reverse in 128 bit groups");
                v = _mm256_permute2x128_si256(v, v, 0x03);                                                         // AVX2 latency = 3
                break; //  low half of first -> hi, high half of second -> low  (coding is 0 1 2 3 => a_lo, a_hi, b_lo, b_hi)  // AVX2
            }

            //        print(r, "r");

            //== now see if we need to shuffle bits.
            if (BIT_GROUP_SIZE < 8) {

              //== next get the lut indices.
              // load constants:
              __m256i _mask_lo = _mm256_load_si256((__m256i*)mask4);                                         // AVX

              //        print(mask_lo, "mask_lo");

              // lower 4 bits
              __m256i lo = _mm256_and_si256(_mask_lo, v); // no and_si256.  had to cast to ps first.           // AVX2
              // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
              // shift not using si128 - that instruction shifts bytes not bits.
              __m256i hi = _mm256_srli_epi16(_mm256_andnot_si256(_mask_lo, v), 4);                                // AVX2

              //        print(lo, "lo");
              //        print(hi, "hi");

              switch (BIT_GROUP_SIZE) {
                case 1:
                  DEBUGF("AVX2 reverse in 1 bit groups");
                  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                  lo = _mm256_shuffle_epi8(_mm256_load_si256((__m256i*)lut1_hi), lo);                        // AVX2
                  hi = _mm256_shuffle_epi8(_mm256_load_si256((__m256i*)lut1_lo), hi);                        // AVX2
                  break;
                case 2:
                  DEBUGF("AVX2 reverse in 2 bit groups");
                  //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
                  lo = _mm256_shuffle_epi8(_mm256_load_si256((__m256i*)lut2_hi), lo);                        // AVX2
                  hi = _mm256_shuffle_epi8(_mm256_load_si256((__m256i*)lut2_lo), hi);                        // AVX2
                  break;
                case 4:
                  DEBUGF("AVX2 reverse in 4 bit groups");
                  lo = _mm256_slli_epi16(lo, 4);  // shift lower 4 to upper										//AVX2
                  break;
              }

              // recombine
              v =  _mm256_or_si256(lo, hi);                                                                  // AVX2

              //      print(res, "result");
            }
            return v;
          }

      };

      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::rev_idx_lane8[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::rev_idx_lane16[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::rev_idx32[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::mask4[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut1_lo[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut1_hi[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut2_lo[32];
      template <unsigned int BIT_GROUP_SIZE, bool POW2>
      constexpr uint8_t bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2, POW2>::lut2_hi[32];

      /// partial template specialization for SSSE3 based bit reverse.  this is defined only for bit_group_sizes that are 1, 2, 4, and 8 (actually powers of 2 up to 128bit)
      template <bool POW2>
      struct bitgroup_ops<3, BIT_REV_AVX2, POW2> {

          /// lookup table for reversing bits in a byte in groups of 4 - no need, just use the mask_lo and shift operator.
          static constexpr unsigned int bitsPerGroup = 3;
          static constexpr unsigned char simd_type = BIT_REV_AVX2;

          static constexpr size_t mask3lo[4] alignas(32) = { 0x9249249249249249, 0x4924924924924924, 0x2492492492492492, 0x9249249249249249 };
          static constexpr size_t mask3mid[4] alignas(32) = { 0x2492492492492492, 0x9249249249249249, 0x4924924924924924, 0x2492492492492492 };
          static constexpr size_t mask3hi[4] alignas(32) = { 0x4924924924924924, 0x2492492492492492, 0x9249249249249249, 0x4924924924924924 };

          bitgroup_ops<1, BIT_REV_AVX2, true> bit_rev_1;

          static ::std::string toString(__m256i const & v) {
            uint8_t tmp[32] alignas(32);
            _mm256_store_si256((__m256i*)tmp, v);
            ::std::stringstream ss;
            for (int i = 31; i >= 0; --i) {
              ss << std::hex << std::setfill('0') << std::setw(2) << static_cast<size_t>(tmp[i]) << " ";
            }
            return ss.str();
          }


          /// shift right by number of bits less than 8.  1 load, 1 shuffle, 2 shifts, 1 or:  5 operations.
          __m256i srli(__m256i val, uint8_t shift) {
            // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
            // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
            //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
            //                        to shift right, b = val, N = 1, lower lane of a should be upper lane of val, upper lane of a should be 0.  possibly 2 ops. permute latency is 3
            __m256i tmp = _mm256_permute2x128_si256(val, val, 0x83);  // lower lane is 0, higher lane is lower lane of val
            tmp = _mm256_alignr_epi8(tmp, val, 1);

            // then left shift to get the bits in the right place
            // shift the input value via epi64 by the number of shifts
            // then or together.
            return _mm256_or_si256(_mm256_slli_epi64(tmp, 8 - shift), _mm256_srli_epi64(val, shift));
          }

          /// shift left by number of bits less than 8.
          __m256i slli(__m256i val, uint8_t shift) {
            // shuffle to get the 8th byte into the 9th position.  shuffle and slli operate on 128 bit lanes only, so do alignr and permute2x128 instead
            // alternative:  alignr.  2 inputs: a, b.  output composition is lowest shift bytes from the 2 lanes of a become the high bytes in the 2 lanes, in order,
            //                              and highest 16-shift bytes from b's 2 lanes are the lower bytes in teh 2 lanes of output.
            //                        to shift left, a = val, N = 15, lower lane of b should be 0, higher lane of b should be lower lane of val.  possibly 2 ops, permute latency is 3.
            __m256i tmp = _mm256_permute2x128_si256(val, val, 0x08);  // lower lane is 0, higher lane is lower lane of val
            tmp = _mm256_alignr_epi8(val, tmp, 15);

            // then left shift to get the bits in the right place
            // shift the input value via epi64 by the number of shifts
            // then or together.
            return _mm256_or_si256(_mm256_srli_epi64(tmp, 8 - shift), _mm256_slli_epi64(val, shift));
          }

          /// reverse function to reverse bits for data types are are not __mm128i  (has to be bigger than at least half)
          size_t reverse(uint8_t * out, uint8_t const * in, size_t len, size_t const bit_offset = 0) {
            if (len == 0) throw ::std::invalid_argument("ERROR reversing byte array: length is 0");

            if (len > 32) throw std::invalid_argument("ERROR:  length of array should be <= 32.  For longer, break it up.");
            //assert(len > 8); //  "ERROR:  length of array should be > 8.  For shorter, use a different bitgroup_ops type.");

            // enforce bit_offset to be less than 8, since we are ORing the first and last bytes.
            if (bit_offset >= 8) throw ::std::invalid_argument("ERROR reversing byte array AVX2:  bit offset should be less than 8.  for larger, shift input instead.");

			// convert to __m128i, then call the __m128i version.  note that there aren't types that have 128 bits other than simd types

            // allocate temp space
            uint8_t tmp[32] alignas(32);
            memset(tmp, 0, 32);

            __m256i v;
            __m256i w;

            // load the data
            if (len == 32) {
            	v = _mm256_loadu_si256((__m256i*)in);
            } else {
            	memcpy(tmp, in, len);

            	v = _mm256_load_si256((__m256i*)tmp);
            }

            // do the conversion
            if (bit_offset == 0) {
            	w = this->reverse( v );
            } else {
            	// ARE THE SHIFTS NECESSARY?  they end up zeroing out the offset, which could be accomplished via a mask.
            	// note that the remainder portion is not zeroed, specifically, if rem == 2, then the center bit is not zeroed, but the side bits are zeroed.
            	// center bit OR with itself does not change.  side bits from next iter ORed with zeroed side bits (this iter) returns the right results.
            	// SO ARE THE ZEROING ALSO NECESSARY?

              // SHIFTS are NOT NECESSARY.



            	// shift right before
//            	v = srli(v, bit_offset);

            	// convert
            	w = this->reverse( v, bit_offset );

            	// shift right after
//            	w = srli(w, bit_offset);
            }

            // copy back to out. need to OR with existing output, so cannot directly store.
            _mm256_store_si256((__m256i*)tmp, w);                    // SSE2

            if (len > 2) memcpy(out+1, tmp + (33 - len), len - 2);
            if (len > 1) out[len - 1] |= tmp[31];
            out[0] |= tmp[32 - len];

            // maskmoveu is very slow.  why?
			//_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));        // SSE2

            // return remainder.
            return (len * 8 - bit_offset) % 3;
          }


          __m256i reverse(__m256i const & u, uint8_t bit_offset = 0 ) {

            // load from memory in reverse is not appropriate here - since we may not have aligned memory, and we have v instead of a memory location.
            // shuffle may be done via register mixing instructions, but may be slower.
            // shuffle of 32 bit and 64 bit may be done via either byte-wise shuffle, or by shuffle_epi32.  should be same speed except that  1 SO post suggests epi32 has lower latency (1 less load due to imm).
            // shuffle of 16 bit can be done via byte-wise shuffle or by 2 calls to shuffle_epi16 plus an or.  probably not faster.
            // CHECK OUT http://www.agner.org/optimize/instruction_tables.pdf for latency and throughput.

            //        print(v, "v");
//            std::cout << "1bit orig: " << toString(u) << std::endl;

            // first reverse bits in groups of 1.
            __m256i v = bit_rev_1.reverse(u);  // 4 loads, 3 shuffles, 3 logical ops, 1 shift op


            // then get the individual bits.  3 loads, 3 ORs.
            __m256i lo = _mm256_and_si256(v, _mm256_load_si256((__m256i*)mask3lo));
            __m256i mid = _mm256_and_si256(v, _mm256_load_si256((__m256i*)mask3mid));
            __m256i hi = _mm256_and_si256(v, _mm256_load_si256((__m256i*)mask3hi));

            __m256i tmp ;

//            std::cout << "bit offset = " << std::dec << static_cast<size_t>(bit_offset) << ", rem = " << static_cast<size_t>((256 - bit_offset) % 3) << std::endl;
            // next shift based on the rem bits.  2 cross boundary shifts = 2 byteshift, 4 shifts, 2 ORs. + 2 ORs.
            if (((256 - bit_offset) % 3) == 2) {
              // rem == 2:                    mid, hi, lo
              tmp  = srli(mid, 2);
              mid = lo;
              lo = tmp;
              hi  = slli( hi, 2);

            } else if (((256 - bit_offset) % 3) == 1) {
              // rem == 1:                    hi, lo, mid
              lo  = srli(lo, 2);
              tmp = slli(mid, 2);
              mid = hi;
              hi  = tmp;
            } else {
              // rem == 0:  first 3 bits are: lo, mid, hi in order of significant bits
              tmp  = srli(hi, 2);
              // mid = mid
              hi  = slli(lo, 2);
              lo = tmp;
            }

            v = _mm256_or_si256(lo, _mm256_or_si256(mid, hi) );

            // alternatively,  for cross-boundary shift, do byte shift, and reverse shift (8-bit_offset), then or. with shifted (offset)->  each CBS is finished in 4 ops w/o load.


            return v;
          }


      };

      template <bool POW2>
      constexpr size_t bitgroup_ops<3, BIT_REV_AVX2, POW2>::mask3lo[4];
      template <bool POW2>
      constexpr size_t bitgroup_ops<3, BIT_REV_AVX2, POW2>::mask3mid[4];
      template <bool POW2>
      constexpr size_t bitgroup_ops<3, BIT_REV_AVX2, POW2>::mask3hi[4];



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

      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is 3, and greater than 0.
       * @param out
       * @param in
       * @param len
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE>
      typename std::enable_if<(BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse(uint8_t * out, uint8_t const * in, size_t const len, uint8_t bit_offset = 0 ) {

        memset(out, 0, len);

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        uint8_t * v = out + len;
        const uint8_t * u = in;


#ifdef __AVX2__
        bitgroup_ops<3, BIT_REV_AVX2, false> op256;
        while (rem >= 32) {
          // enough bytes.  do an iteration
          v -= 32;
          bit_offset = op256.reverse(v, u, 32, bit_offset);
          u += 32;
          rem -= 32;

          if (bit_offset > 0) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
//          std::cout << "AVX full.  rem=" << std::dec << rem << " len=" << (len - rem) << std::endl;
        }
        if (rem > 16) {
          // do another iteration with all the remaining.   not that offset is tied to the current u position, so can't use memcpy avoiding approach.
          bit_offset = op256.reverse(out, in + len - rem, rem, bit_offset);
          rem = 0;
//          std::cout << "AVX part.  rem=" << std::dec << rem << " len=" << (len - rem) << std::endl;
        }  // otherwise use either SSSE3 or SWAR
#endif

#ifdef __SSSE3__
        bitgroup_ops<3, BIT_REV_SSSE3, false> op128;
        while (rem >= 16) {

          // enough bytes.  do an iteration
          v -= 16;
          bit_offset = op128.reverse(v, u, 16, bit_offset);
          u += 16;
          rem -= 16;

          if (bit_offset > 0) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
//          std::cout << "SSSE full.  rem = " << std::dec <<  rem << " len=" << (len - rem) << std::endl;
        }
        if (rem > 8) {

          // do another iteration with all the remaining.
          bit_offset = op128.reverse(out, in + len - rem, rem, bit_offset);
          rem = 0;
//          std::cout << "SSSE part.  rem = " << std::dec << rem << " len=" << (len - rem) << std::endl;
        }
#endif
        bitgroup_ops<3, BIT_REV_SWAR> op64;
        while (rem >= 8) {

          // enough bytes.  do an iteration
          v -= 8;
          bit_offset = op64.reverse(v, u, 8, bit_offset);
          u += 8;
          rem -= 8;

          if (bit_offset > 0) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
//          std::cout << "SWAR full.  rem = " << std::dec << rem << " len=" << (len - rem) << std::endl;

        }
        if (rem > 0) {

          // do another iteration with all the remaining.

          // original length is less than 32, so has to do memcpy
          bit_offset = op64.reverse(out, in + len - rem, rem, bit_offset);

          rem = 0;
//          std::cout << "SWAR part.  rem = " << std::dec << rem << " len=" << (len - rem) << std::endl;
        }

        return bit_offset;  // return remainder.
      }

      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE>
      typename std::enable_if<((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse(uint8_t * out, uint8_t const * in, size_t const len, uint8_t bit_offset = 0 ) {
        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= 64, "ERROR: currenly reverse does not support 128 BIT_GRUOP_SIZE");

        if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE require bit_offset to be 0.");

        if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");

        memset(out, 0, len);

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        uint8_t * v = out + len;
        const uint8_t * u = in;

#ifdef __AVX2__
        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_AVX2> op256;
        while (rem >= 32) {
          // enough bytes.  do an iteration
          v -= 32;
          op256.reverse(v, u, 32, 0);
          u += 32;

          rem -= 32;
        }
        if (rem > 16) {
          // do another iteration with all the remaining.
          if (len > 32) {  // original length has 32 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            op256.reverse(out, in + len - 32, 32, 0);
          } else {  // original length is less than 32, so has to do memcpy
            op256.reverse(out, in + len - rem, rem, 0);
          }
          rem = 0;
        }  // otherwise use either SSSE3 or SWAR
#endif

#ifdef __SSSE3__
        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SSSE3> op128;
        while (rem >= 16) {
          // enough bytes.  do an iteration
          v -= 16;
          op128.reverse(v, u, 16, 0);
          u += 16;

          rem -= 16;
        }
        if (rem > 8) {
          // do another iteration with all the remaining.
          if (len > 16) {  // original length has 32 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            op128.reverse(out, in + len - 16, 16, 0);
          } else {  // original length is less than 32, so has to do memcpy
            op128.reverse(out, in + len - rem, rem, 0);
          }
          rem = 0;
        }
#endif
        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SWAR> op64;
        while (rem >= 8) {
          // enough bytes.  do an iteration
          v -= 8;
          op64.reverse(v, u, 8, 0);
          u += 8;

          rem -= 8;
        }
        if (rem > 0) {
          // do another iteration with all the remaining.
          if (len > 8) {  // original length has 32 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            op64.reverse(out, in + len - 8, 8, 0);
          } else {  // original length is less than 32, so has to do memcpy
            op64.reverse(out, in + len - rem, rem, 0);
          }
          rem = 0;
        }

        return 0;  // return remainder.
      }


      /// convenience function to convert from WORD_TYPE pointer to uint8_t pointer.
      template <typename WORD_TYPE, unsigned int BIT_GROUP_SIZE, typename ::std::enable_if<(::std::is_integral<WORD_TYPE>::value && (sizeof(WORD_TYPE) > 1)), int>::type = 1>
      uint8_t reverse(WORD_TYPE * out, WORD_TYPE const * in, size_t const len, uint8_t const bit_offset) {
        return reverse<BIT_GROUP_SIZE>(reinterpret_cast<uint8_t*>(out), reinterpret_cast<uint8_t const *>(in), len * sizeof(WORD_TYPE), bit_offset);
      }


      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is 3, and greater than 0.
       * @param out
       * @param in
       * @param len
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE>
      typename std::enable_if<(BIT_GROUP_SIZE == 3), uint8_t>::type
      reverse_seq(uint8_t * out, uint8_t const * in, size_t const len, uint8_t bit_offset = 0 ) {

        memset(out, 0, len);

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        uint8_t * v = out + len;
        const uint8_t * u = in;


        bitgroup_ops<3, BIT_REV_SEQ> op64;
        while (rem >= 8) {

          // enough bytes.  do an iteration
          v -= 8;
          bit_offset = op64.reverse(v, u, 8, bit_offset);
          u += 8;
          rem -= 8;

          if (bit_offset > 0) {  // if there is overlap.  adjust
            u -= 1;
            rem += 1;
            v += 1;
            bit_offset = 8 - bit_offset;
          } // else no adjustment is needed.
//          std::cout << "SWAR full.  rem = " << std::dec << rem << " len=" << (len - rem) << std::endl;

        }
        if (rem > 0) {

          // do another iteration with all the remaining.

          // original length is less than 32, so has to do memcpy
          bit_offset = op64.reverse(out, in + len - rem, rem, bit_offset);

          rem = 0;
//          std::cout << "SWAR part.  rem = " << std::dec << rem << " len=" << (len - rem) << std::endl;
        }

        return bit_offset;  // return remainder.
      }

      /**
       * @brief
       * @details   enabled only if BIT_GROUP_SIZE is power of 2, and greater than 0.
       * @param out
       * @param in
       * @param len
       * @param bit_offset
       * @return
       */
      template <unsigned int BIT_GROUP_SIZE>
      typename std::enable_if<((BIT_GROUP_SIZE & (BIT_GROUP_SIZE - 1)) == 0), uint8_t>::type
      reverse_seq(uint8_t * out, uint8_t const * in, size_t const len, uint8_t bit_offset = 0 ) {
        static_assert(BIT_GROUP_SIZE > 0, "ERROR: BIT_GROUP_SIZE cannot be 0");
        static_assert(BIT_GROUP_SIZE <= 64, "ERROR: currenly reverse does not support 128 BIT_GRUOP_SIZE");

        if (bit_offset != 0) throw std::invalid_argument("ERROR: power of 2 BIT_GROUP_SIZE require bit_offset to be 0.");

        if ((len % ((BIT_GROUP_SIZE + 7) / 8)) > 0)
          throw ::std::invalid_argument("ERROR reversing byte array:  if BIT_GROUP_SIZE > 8 bits, len needs to be a multiple of BIT_GROUP_SIZE in bytes");

        memset(out, 0, len);

        // BIT_GROUP_SIZE is an accelerated one.  so decide based on len and available instruction sets.
        size_t rem = len;
        // pointers
        uint8_t * v = out + len;
        const uint8_t * u = in;

        bitgroup_ops<BIT_GROUP_SIZE, BIT_REV_SEQ> op64;
        while (rem >= 8) {
          // enough bytes.  do an iteration
          v -= 8;
          op64.reverse(v, u, 8, 0);
          u += 8;

          rem -= 8;
        }
        if (rem > 0) {
          // do another iteration with all the remaining.
          if (len > 8) {  // original length has 32 bytes or more, so avoid memcpy.  duplicate a little work but that's okay.
            op64.reverse(out, in + len - 8, 8, 0);
          } else {  // original length is less than 32, so has to do memcpy
            op64.reverse(out, in + len - rem, rem, 0);
          }
          rem = 0;
        }

        return 0;  // return remainder.
      }


      /// convenience function to convert from WORD_TYPE pointer to uint8_t pointer.
      template <typename WORD_TYPE, unsigned int BIT_GROUP_SIZE, typename ::std::enable_if<(::std::is_integral<WORD_TYPE>::value && (sizeof(WORD_TYPE) > 1)), int>::type = 1>
      uint8_t reverse_seq(WORD_TYPE * out, WORD_TYPE const * in, size_t const len, uint8_t const bit_offset) {
        return reverse_seq<BIT_GROUP_SIZE>(reinterpret_cast<uint8_t*>(out), reinterpret_cast<uint8_t const *>(in), len * sizeof(WORD_TYPE), bit_offset);
      }



    } // namespace bit_ops


  } // namespace utils

} // namespace bliss








#endif /* SRC_UTILS_BITGROUP_OPS_HPP_ */
