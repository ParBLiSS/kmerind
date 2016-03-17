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
 * common.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: tpan
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <random>
#include <array>
#include <cstdint>
#include <utility>
#include <iostream>
#include <tuple>
#include <limits>

// include files to test
#include "utils/bitgroup_ops.hpp"
#include "utils/test/bit_test_common.hpp"

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bit reverse.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.


// NOTE if the gtest fixture class is missing the trailing semicolon, compile error about "expected initializer..."
template <unsigned char BITS_PER_GROUP>
class BitReverseTestHelper {
  public:

    uint8_t input[32] = {  0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10 };

    uint8_t array[128] = { 0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10};


    // bit group size larger than byte
    template <unsigned char BITS = BITS_PER_GROUP, typename std::enable_if<(BITS >= 8), int>::type = 1>
    bool is_reverse(uint8_t const * orig, uint8_t const* rev, size_t len, uint8_t bit_offset = 0) {
      bool same = true;
      bool local_same = false;

      //std::cout << "u: " << std::hex << orig << " v: " << std::hex << rev << std::endl;


      // compare based on BITS size.
      // reverse bytes or larger, than just compare the bytes.
      const uint8_t* w = rev + len;
      const uint8_t* uu = orig;

      uint8_t bytesPerChar = BITS / 8;
      uint8_t max = len / bytesPerChar;

      for (size_t i = 0; i < max; ++i) {
        w -= bytesPerChar;
        local_same = (memcmp(uu, w, bytesPerChar) == 0);

        if (!local_same) {
          std::cout << "len " << std::dec << len << " not matched at " << static_cast<size_t>(i) << " w: ";
          for (int k = 0; k < bytesPerChar; ++k) {
            std::cout << std::hex << static_cast<size_t>(w[bytesPerChar - 1 - k]) << " ";
          }
          std::cout << "uu: ";
          for (int k = 0; k < bytesPerChar; ++k) {
            std::cout << std::hex << static_cast<size_t>(uu[bytesPerChar - 1 - k]) << " ";
          }
          std::cout << std::endl;
        }

        same &= local_same;
        uu += bytesPerChar;
      }

      return same;
    }

    // power of 2 bit group size, bit group smaller than byte
    template <unsigned char BITS = BITS_PER_GROUP, typename std::enable_if<(BITS < 8) && ((BITS & (BITS - 1)) == 0), int>::type = 1>
    bool is_reverse(uint8_t const * orig, uint8_t const* rev, size_t len, uint8_t bit_offset = 0) {
      bool same = true;


      const uint8_t* uu = orig;
      const uint8_t* w = rev + len;

      uint8_t bits_mask = static_cast<uint8_t>(~(0xFF << BITS));

      for (size_t i = 0; i < len; ++i) {
        --w;

        for (uint8_t j = 0, k = 8 - BITS; j < 8; j += BITS, k -= BITS) {
          same &= ((*uu >> j) & bits_mask) == ((*w >> k) & bits_mask);
        }

        ++uu;
      }

      return same;
    }

    // non power of 2 bit group.  bit group size less than 8
    template <unsigned char BITS = BITS_PER_GROUP, typename std::enable_if<(BITS < 8) && ((BITS & (BITS - 1)) != 0), int>::type = 1>
    bool is_reverse(uint8_t const * orig, uint8_t const* rev, size_t len, uint8_t bit_offset = 0) {
      bool same = true;
      bool local_same = false;

      if (len == 1) {  // single byte

        size_t rem = (8 - bit_offset) % BITS;
        unsigned char bits_mask = static_cast<unsigned char>(~(0xFF << BITS));
        //unsigned char rem_mask = static_cast<unsigned char>(~(0xFF >> rem));

        // check the last part is same.
        //same = (*orig & rem_mask) == (*rev & rem_mask);

        // check the reversed part
        // align test data so the reversed part is at the highest bits
        for (size_t i = bit_offset, j = (8 - bit_offset - BITS); i < (8 - rem); i += BITS, j -= BITS) {
          local_same = ((*orig >> i) & bits_mask) == ((*rev >> j) & bits_mask);

          same &= local_same;
        }

      } else {   // multibyte

        // use 2 bytes, overlapped by 1 byte in successive run, then all bits will be covered.

        uint8_t offset = bit_offset;  // offset from lowest.

        const uint8_t* w = rev + len - 1;  // doing it in pairs of 2.
        const uint8_t* uu = orig;

        uint16_t bits_mask = static_cast<uint16_t>(~(0xFFFF << BITS));


        for (uint8_t i = 0; i < len-1; ++i) {
          --w;

          size_t rem = (16 - offset) % BITS;  // remainder within the current short.

          for (size_t k = 16 - offset - BITS; offset < (16 - rem); offset += BITS, k -= BITS) {
            local_same = ((*(reinterpret_cast<const uint16_t*>(uu)) >> offset) & bits_mask) == ((*(reinterpret_cast<const uint16_t*>(w)) >> k) & bits_mask);

            if (!local_same) {
              std::cout << "diff @ " << std::dec << static_cast<size_t>(i) << ". u_off=" << static_cast<size_t>(offset) << ", v_offset=" << static_cast<size_t>(k) << ":";
              std::cout << " full u=" << std::hex << std::setw(2 * sizeof(uint16_t)) << std::setfill('0') << *(reinterpret_cast<const uint16_t*>(uu)) << " full v=" << std::setw(2 * sizeof(uint16_t)) << std::setfill('0') << *(reinterpret_cast<const uint16_t*>(w));
              std::cout << " u=" << std::hex << std::setw(2 * sizeof(uint16_t)) << std::setfill('0') << static_cast<size_t>(((*(reinterpret_cast<const uint16_t*>(uu)) >> offset) & bits_mask));
              std::cout << " v=" << std::hex << std::setw(2 * sizeof(uint16_t)) << std::setfill('0') << static_cast<size_t>(((*(reinterpret_cast<const uint16_t*>(w)) >> k) & bits_mask)) << std::endl;
            }
            same &= local_same;
          }
          // note:  offset is < 16-rem inside the loop, at 16 - rem when exiting.  new offset is shifted by 8.  not same as modulus. e.g. when rem is 0, the next offset should be 8.
          offset = (16-rem) - 8;

          ++uu;
        }


      }

      return same;

    }

    // enable_if is necessary.  else WORD_TYPE may match to pointer types.
    template <typename WORD_TYPE, typename std::enable_if<std::is_integral<WORD_TYPE>::value, int>::type = 1>
    bool is_reverse(WORD_TYPE const & orig, WORD_TYPE const & rev, uint8_t bit_offset = 0) {
      return is_reverse(reinterpret_cast<const uint8_t*>(&orig),
                        reinterpret_cast<const uint8_t*>(&rev), sizeof(WORD_TYPE), bit_offset);
    }

#ifdef __SSSE3__
    bool is_reverse(__m128i const & orig, __m128i const & rev, uint8_t bit_offset = 0) {

      uint8_t BLISS_ALIGNED_ARRAY(orig_arr, 16, 16);
      uint8_t  BLISS_ALIGNED_ARRAY(rev_arr, 16, 16);

      _mm_store_si128((__m128i*)orig_arr, orig);
      _mm_store_si128((__m128i*)rev_arr, rev);

      return is_reverse(orig_arr,
                        rev_arr, 16, bit_offset);
    }
#endif

#ifdef __AVX2__

    bool is_reverse(__m256i const & orig, __m256i const & rev, uint8_t bit_offset = 0) {

      uint8_t BLISS_ALIGNED_ARRAY(orig_arr, 32, 32);
      uint8_t  BLISS_ALIGNED_ARRAY(rev_arr, 32, 32);

      _mm256_store_si256((__m256i*)orig_arr, orig);
      _mm256_store_si256((__m256i*)rev_arr, rev);

      return is_reverse(orig_arr,
                        rev_arr, 32, bit_offset);
    }
#endif
};





#endif /* COMMON_HPP_ */
