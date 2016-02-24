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


// include google test
#include <gtest/gtest.h>

#include <random>
#include <array>
#include <cstdint>
#include <utility>
#include <iostream>
#include <tuple>
#include <limits>

// include files to test
#include "utils/bitgroup_ops.hpp"

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bit reverse.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.

template <unsigned char Bits>
struct BitsParam { static constexpr unsigned char bitsPerGroup = Bits; };



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

      for (uint8_t i = 0; i < max; ++i) {
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

      for (uint8_t i = 0; i < len; ++i) {
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

          for (uint8_t k = 16 - offset - BITS; offset < (16 - rem); offset += BITS, k -= BITS) {
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


template <typename T>
class BitReverseTest : public ::testing::Test {
  protected:

    BitReverseTestHelper<T::bitsPerGroup> helper;

    template <typename WORD_TYPE>
    bool is_reverse(WORD_TYPE const & orig, WORD_TYPE const & rev, uint8_t bit_offset = 0) {

      return helper.is_reverse(orig, rev, bit_offset);
    }

    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len, uint8_t bit_offset = 0)  {

      return helper.is_reverse(out, in, len, bit_offset);
    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseTest);

TYPED_TEST_P(BitReverseTest, reverse_uint8)
{
  using WORD_TYPE = uint8_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

      WORD_TYPE in;
  	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));


      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_uint16)
{
  using WORD_TYPE = uint16_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

        WORD_TYPE in;
    	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));

      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_uint32)
{
  using WORD_TYPE = uint32_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

        WORD_TYPE in;
    	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));

      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}


TYPED_TEST_P(BitReverseTest, reverse_uint64)
{
  using WORD_TYPE = uint64_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

        WORD_TYPE in;
    	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));

      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_short_array)
{

  TypeParam op;

  uint8_t BLISS_ALIGNED_ARRAY(out, 32, 32);

  unsigned int max = 8;

  for (unsigned int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


    if (TypeParam::bitsPerGroup == 3) {
      for (unsigned int k = 0; k <= (32 - i); ++k) {
        for (uint8_t j = 0; j < 8; ++j) {
          memset(out, 0, 32);

          op.reverse(out, this->helper.input + k, i, j);

          bool same = this->is_reverse(this->helper.input + k, out, i, j);

          if (!same) {

            printf("array size = %d\n", i);
          }

          EXPECT_TRUE(same);
        }
      }
    } else {
      for (unsigned int k = 0; k <= (32 - i); ++k) {
        memset(out, 0, 32);

        op.reverse(out, this->helper.input + k, i);

        bool same = this->is_reverse(this->helper.input + k, out, i);

        if (!same) {
          std::cout << "in (MSB to LSB): ";
          for (int64_t j = i-1; j >= 0; --j) {
            std::cout << std::hex << static_cast<size_t>(this->helper.input[k + i - 1 - j]) << " ";

          }
          std::cout << std::endl;

          std::cout << "out (MSB to LSB): ";
          for (int64_t j = i-1; j >= 0; --j) {
            std::cout << std::hex << static_cast<size_t>(out[i - 1 - j]) << " ";

          }
          std::cout << std::endl;

          printf("array size = %d\n", i);
        }

        EXPECT_TRUE(same);
      }
    }
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseTest,
                           reverse_uint8,
                           reverse_uint16,
                           reverse_uint32,
                           reverse_uint64,
                           reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
     ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
      ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
       ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
        ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
         ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
          ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
           ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
            ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
             ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
              ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
               ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
                ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
                 ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_SWAR>
> BitReverseTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseTest, BitReverseTestTypes);




#ifdef __SSSE3__

template <typename T>
class BitReverseSSSETest : public ::testing::Test {
  protected:

    BitReverseTestHelper<T::bitsPerGroup> helper;

    bool is_reverse(__m128i const & orig, __m128i const & rev, uint8_t bit_offset = 0) {

      return helper.is_reverse(orig, rev, bit_offset);
    }

    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len, uint8_t bit_offset = 0)  {

      return helper.is_reverse(out, in, len, bit_offset);
    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseSSSETest);

TYPED_TEST_P(BitReverseSSSETest, reverse_m128i)
{
  TypeParam op;

  if (TypeParam::bitsPerGroup < 128 ) {
    for (int k = 0; k < 16; ++k) {
      __m128i in = _mm_loadu_si128((__m128i*)(this->helper.input + k));

      __m128i out = op.reverse(in);

      EXPECT_TRUE(this->is_reverse(in, out, 0));
    }
  }  // else too large, so don't do the test.

}


TYPED_TEST_P(BitReverseSSSETest, reverse_short_array)
{

  TypeParam op;

  uint8_t BLISS_ALIGNED_ARRAY(out, 32, 32);

  int max = 16;

  for (int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


    if (TypeParam::bitsPerGroup == 3) {
      for (int k = 0; k <= (32 - i); ++k) {
        for (int j = 0; j < 8; ++j) {


          memset(out, 0, 32);

          op.reverse(out, this->helper.input + k, i, j);

          bool same = this->is_reverse(this->helper.input + k, out, i, j);

          if (!same) {

            printf("array size = %d\n", i);
          }

          EXPECT_TRUE(same);
        }
      }
    } else {

      for (int k = 0; k <= (32 - i); ++k) {
        memset(out, 0, 32);

        op.reverse(out, this->helper.input + k, i);

        bool same = this->is_reverse(this->helper.input + k, out, i);

        if (!same) {

          printf("array size = %d\n", i);
        }

        EXPECT_TRUE(same);
      }


    }
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseSSSETest, reverse_m128i, reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
     ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
      ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
       ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
        ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
         ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
          ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
           ::bliss::utils::bit_ops::bitgroup_ops<64, ::bliss::utils::bit_ops::BIT_REV_SSSE3>
> BitReverseSSSETestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseSSSETest, BitReverseSSSETestTypes);

#endif



#ifdef __AVX2__
template <typename T>
class BitReverseAVX2Test : public ::testing::Test {
  protected:

    BitReverseTestHelper<T::bitsPerGroup> helper;

    bool is_reverse(__m256i const & orig, __m256i const & rev, uint8_t bit_offset = 0) {

      return helper.is_reverse(orig, rev, bit_offset);
    }

    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len, uint8_t bit_offset = 0)  {

      return helper.is_reverse(out, in, len, bit_offset);
    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseAVX2Test);

TYPED_TEST_P(BitReverseAVX2Test, reverse_m256i)
{
  TypeParam op;

  if (TypeParam::bitsPerGroup < 256 ) {
    __m256i in = _mm256_loadu_si256((__m256i*)(this->helper.input));

    __m256i out = op.reverse(in);

    ASSERT_TRUE(this->is_reverse(in, out));
  }  // else too large, so don't do the test.

}


TYPED_TEST_P(BitReverseAVX2Test, reverse_short_array)
{

  TypeParam op;

  uint8_t BLISS_ALIGNED_ARRAY(out, 32, 32);

  int max = 32;


  for (int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


    if (TypeParam::bitsPerGroup == 3) {
      for (int k = 0; k <= (32 - i); ++k) {
        for (int j = 0; j < 8; ++j) {
          memset(out, 0, 32);

          op.reverse(out, this->helper.input + k, i, j);

          bool same = this->is_reverse(this->helper.input + k, out, i, j);

          if (!same) {
            std::cout << "in: ";
            for (int l = 0; l < i; ++l) {
              std::cout << std::hex << static_cast<size_t>(this->helper.input[k + i - 1 - l]) << " ";

            }
            std::cout << std::endl;

            std::cout << "out: ";
            for (int l = 0; l < i; ++l) {
              std::cout << std::hex << static_cast<size_t>(out[i - 1 - l]) << " ";

            }
            std::cout << std::endl;


            printf("array size = %d, offset = %d, input byte offset = %d \n", i, j, k);
          }

          ASSERT_TRUE(same);
        }
      }
    } else {
      for (int k = 0; k <= (32 - i); ++k) {
        memset(out, 0, 32);

        op.reverse(out, this->helper.input + k, i, 0);

        bool same = this->is_reverse(this->helper.input + k, out, i, 0);

        if (!same) {
          printf("array size = %d\n", i);
        }

        ASSERT_TRUE(same);
      }
    }
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseAVX2Test, reverse_m256i, reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
     ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
      ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
       ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
        ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
         ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
          ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
           ::bliss::utils::bit_ops::bitgroup_ops<64, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
            ::bliss::utils::bit_ops::bitgroup_ops<128, ::bliss::utils::bit_ops::BIT_REV_AVX2>
> BitReverseAVX2TestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseAVX2Test, BitReverseAVX2TestTypes);

#endif




//=========== test long array


template <typename T>
class BitReverseLongArrayTest : public ::testing::Test {
  protected:

    BitReverseTestHelper<T::bitsPerGroup> helper;

    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len, uint8_t bit_offset = 0)  {

      return helper.is_reverse(out, in, len, bit_offset);
    }

};

TYPED_TEST_CASE_P(BitReverseLongArrayTest);

TYPED_TEST_P(BitReverseLongArrayTest, reverse_long_array)
{

  uint8_t BLISS_ALIGNED_ARRAY(out, 128, 32);

  unsigned int max = 128;

  for (unsigned int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


      for (unsigned int k = 0; k <= (128 - i); ++k) {
        memset(out, 0, 128);

        bliss::utils::bit_ops::reverse<TypeParam::bitsPerGroup, bliss::utils::bit_ops::BIT_REV_AVX2>(out, this->helper.input + k, i);

        bool same = this->is_reverse(this->helper.input + k, out, i);

        if (!same) {

          printf("array size = %d\n", i);
        }

        EXPECT_TRUE(same);
      }
  }
}

REGISTER_TYPED_TEST_CASE_P(BitReverseLongArrayTest, reverse_long_array);



typedef ::testing::Types<
    BitsParam< 1> ,
    BitsParam< 2> ,
    BitsParam< 3> ,
    BitsParam< 4> ,
    BitsParam< 8> ,
    BitsParam<16> ,
    BitsParam<32>
> BitReverseLongArrayTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseLongArrayTest, BitReverseLongArrayTestTypes);





template <typename T>
class BitReverseFixedArrayTest : public ::testing::Test {
  protected:

    template <uint8_t BIT_GROUP_SIZE>
    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len)  {
      BitReverseTestHelper<BIT_GROUP_SIZE> helper;

      return helper.is_reverse(out, in, len, 0);
//                               (len * 8) % BIT_GROUP_SIZE);  // if BITS = 3, offset needs to be calculated.
    }

    template <unsigned int BITS, uint8_t MAX_SIMD_TYPE,
      typename data_type,
      size_t data_size,
      typename ::std::enable_if<((data_size * sizeof(data_type) * 8) <= BITS), int>::type = 1>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {
      if ((data_size * sizeof(data_type) * 8) == BITS) {
        memcpy(out, in, data_size * sizeof(data_type)); return;
      } else if ((data_size * sizeof(data_type) * 8) < BITS) {
        return;
      }
    }
    template <unsigned int BITS, uint8_t MAX_SIMD_TYPE,
      typename data_type,
      size_t data_size,
      typename ::std::enable_if<((data_size * sizeof(data_type) * 8) > BITS), int>::type = 1>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {

      BitReverseTestHelper<BITS> helper;
      bool same = true;

      for (unsigned int k = 0; k < (32 - data_size); ++k) {
//          printf("array size = %lu, sizeof(datatype) = %lu, bits = %u, SIMD = %u, k = %u\n", data_size, sizeof(data_type), BITS, MAX_SIMD_TYPE, k);

        memcpy(in, helper.array + k, data_size * sizeof(data_type));

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template reverse<BITS, MAX_SIMD_TYPE>(out, in);

        same =
            this->is_reverse<BITS>(helper.array + k,
                                   reinterpret_cast<uint8_t*>(out),
                                   data_size * sizeof(data_type));

        if (!same) {
          std::cout << "in (MSB to LSB): ";
          for (int64_t j = data_size - 1; j >= 0; --j) {
            std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(in[j]) << " ";
          }
          std::cout << std::endl;

          std::cout << "out (MSB to LSB): ";
          for (int64_t j = data_size - 1; j >= 0; --j) {
            std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(out[j]) << " ";

          }
          std::cout << std::endl;

        }

        ASSERT_TRUE(same);

      }
    }

    template <uint8_t SIMDType, typename data_type, size_t data_size>
    void run_tests() {
    	  data_type in  [data_size];
    	  data_type out [data_size];
    	  // NEED TO SPECIFY data_type and data_size for icc.  not needed for clang or gcc.
    	  switch (sizeof(data_type)) {
    	    case 8:
    	      test<32, SIMDType, data_type, data_size>(in, out);
    	    case 4:
    	      test<16, SIMDType, data_type, data_size>(in, out);
    	    case 2:
    	      test<8,  SIMDType, data_type, data_size>(in, out);
    	    case 1:
    	      test<4,  SIMDType, data_type, data_size>(in, out);
    	      test<2,  SIMDType, data_type, data_size>(in, out);
    	      test<1,  SIMDType, data_type, data_size>(in, out);
    	      test<3,  SIMDType, data_type, data_size>(in, out);
    	      break;
    	    default:
    	      break;
    	  }

    }
};

TYPED_TEST_CASE_P(BitReverseFixedArrayTest);

TYPED_TEST_P(BitReverseFixedArrayTest, reverse_seq)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_SEQ, data_type, data_size>();
}


TYPED_TEST_P(BitReverseFixedArrayTest, reverse_swar)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_SWAR, data_type, data_size>();
}

#ifdef __SSSE3__
TYPED_TEST_P(BitReverseFixedArrayTest, reverse_ssse3)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_SSSE3, data_type, data_size>();

}
#endif

#ifdef __AVX2__
TYPED_TEST_P(BitReverseFixedArrayTest, reverse_avx2)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_AVX2, data_type, data_size>();
}
#endif


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseFixedArrayTest,
#ifdef __SSSE3__
                           reverse_ssse3,
#endif
#ifdef __AVX2__
                           reverse_avx2,
#endif
                           reverse_swar,
                           reverse_seq);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
	    ::std::tuple<BitsParam< 1>, uint8_t >,    //  1  (uint8_t)
		::std::tuple<BitsParam< 1>, uint16_t >,   //  2  (uint16_t)
		::std::tuple<BitsParam< 3>, uint8_t >,    //  3  1.5 uint16_t
		::std::tuple<BitsParam< 1>, uint32_t >,   //  4  (uint32_t)
		::std::tuple<BitsParam< 3>, uint16_t >,   //  6  1.5 uint32_t
		::std::tuple<BitsParam< 1>, uint64_t >,   //  8   (uint64_t)
		::std::tuple<BitsParam< 3>, uint32_t >,   //  12  1.5 uint64_t
		::std::tuple<BitsParam< 2>, uint64_t >,   //  16  (__m128i)
		::std::tuple<BitsParam< 3>, uint64_t >,   //  24  1.5 __m128i
		::std::tuple<BitsParam< 4>, uint64_t >,   //  32  (__m256i)
		::std::tuple<BitsParam< 5>, uint64_t >,   //  40  1.5 __m256i
		::std::tuple<BitsParam< 8>, uint64_t >    //  64  (512 bits)
> BitReverseFixedArrayTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseFixedArrayTest, BitReverseFixedArrayTestTypes);









template <typename T>
class BitReverseTransformTest : public ::testing::Test {
  protected:

    template <uint8_t BIT_GROUP_SIZE,  uint8_t pad_bits,
    	typename data_type,
    	size_t data_size >
    bool is_reverse(data_type (&in)[data_size], data_type (&out)[data_size])  {
      BitReverseTestHelper<BIT_GROUP_SIZE> helper;

      // shift input left by padbits
      using SIMD = bliss::utils::bit_ops::BITREV_AUTO_AGGRESSIVE<(sizeof(data_type) * data_size), bliss::utils::bit_ops::BITREV_SWAR>;

      data_type tmp[data_size];
      bliss::utils::bit_ops::left_shift<SIMD, pad_bits>(tmp, out);

      // clear the pad bits in input

      data_type tmp2[data_size];
      memcpy(tmp2, in, sizeof(data_type) * data_size);
      tmp2[data_size - 1] &= (::std::numeric_limits<data_type>::max() >> pad_bits);

      // then compare normally.
      return helper.is_reverse(reinterpret_cast<uint8_t const *>(tmp), reinterpret_cast<uint8_t const *>(tmp2), (sizeof(data_type) * data_size));

    }

    template <unsigned int BITS, typename MAX_SIMD_TYPE, uint8_t pad_bits,
      typename data_type,
      size_t data_size,
      typename ::std::enable_if<((data_size * sizeof(data_type) * 8) % BITS > 0) || (pad_bits >= (sizeof(data_type) * 8)), int>::type = 1>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {
      if ((data_size * sizeof(data_type) * 8) == BITS) {
        memcpy(out, in, data_size * sizeof(data_type)); return;
      } else if ((data_size * sizeof(data_type) * 8) < BITS) {
        return;
      }
    }
    template <unsigned int BITS, typename MAX_SIMD_TYPE, uint8_t pad_bits,
      typename data_type,
      size_t data_size,
      typename ::std::enable_if<((data_size * sizeof(data_type) * 8) % BITS == 0) && (pad_bits < (sizeof(data_type) * 8)), int>::type = 1>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {

      BitReverseTestHelper<BITS> helper;

      for (unsigned int k = 0; k < (32 - data_size); ++k) {

        memcpy(in, helper.array + k, data_size * sizeof(data_type));

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template reverse<BITS, MAX_SIMD_TYPE, pad_bits>(out, in);

        bool same = this->is_reverse<BITS, pad_bits>(in, out);

        if (!same) {
          std::cout << "in (MSB to LSB): ";
          for (int64_t j = data_size - 1; j >= 0; --j) {
            std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(in[j]) << " ";
          }
          std::cout << std::endl;

          std::cout << "out (MSB to LSB): ";
          for (int64_t j = data_size - 1; j >= 0; --j) {
            std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(out[j]) << " ";

          }
          std::cout << std::endl;

          printf("array size = %lu, sizeof(datatype) = %lu, bits = %u, SIMD = %u, k = %u\n", data_size, sizeof(data_type), BITS, MAX_SIMD_TYPE::SIMDVal, k);
        }

        ASSERT_TRUE(same);


      }
    }

    template <typename SIMDType, typename data_type, size_t data_size>
    void run_tests() {
    	  data_type in  [data_size];
    	  data_type out [data_size];
    	  // NEED TO SPECIFY data_type and data_size for icc.  not needed for clang or gcc.
    	  switch (sizeof(data_type)) {
    	    case 8:
    	      test<32, SIMDType,  0>(in, out);
    	      test<32, SIMDType, 32>(in, out);
    	    case 4:
    	      test<16, SIMDType,  0>(in, out);
    	      test<16, SIMDType, 16>(in, out);
    	    case 2:
    	      test<8,  SIMDType, 0>(in, out);
    	      test<8,  SIMDType, 8>(in, out);
    	    case 1:
    	      test<4,  SIMDType,  0>(in, out);
    	      test<4,  SIMDType,  4>(in, out);
    	      test<2,  SIMDType,  0>(in, out);
    	      test< 2, SIMDType,  2>(in, out);
    	      test< 2, SIMDType,  6>(in, out);
    	      test<1,  SIMDType,  0>(in, out);
    	      test< 1, SIMDType,  1>(in, out);
    	      test< 1, SIMDType,  7>(in, out);
    	      test< 3, SIMDType,  ((sizeof(data_type) * data_size * 8) % 3 )>(in, out);
    	      break;
    	    default:
    	      break;
    	  }

    }
};

TYPED_TEST_CASE_P(BitReverseTransformTest);

TYPED_TEST_P(BitReverseTransformTest, reverse_seq)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SEQ, data_type, data_size>();
}


TYPED_TEST_P(BitReverseTransformTest, reverse_swar)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size>();
}

#ifdef __SSSE3__
TYPED_TEST_P(BitReverseTransformTest, reverse_ssse3)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size>();

}
#endif

#ifdef __AVX2__
TYPED_TEST_P(BitReverseTransformTest, reverse_avx2)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size>();
}
#endif


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseTransformTest,
#ifdef __SSSE3__
                           reverse_ssse3,
#endif
#ifdef __AVX2__
                           reverse_avx2,
#endif
                           reverse_swar,
                           reverse_seq);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
	    ::std::tuple<BitsParam< 1>, uint8_t >,    //  1  (uint8_t)
		::std::tuple<BitsParam< 1>, uint16_t >,   //  2  (uint16_t)
		::std::tuple<BitsParam< 3>, uint8_t >,    //  3  1.5 uint16_t
		::std::tuple<BitsParam< 1>, uint32_t >,   //  4  (uint32_t)
		::std::tuple<BitsParam< 3>, uint16_t >,   //  6  1.5 uint32_t
		::std::tuple<BitsParam< 1>, uint64_t >,   //  8   (uint64_t)
		::std::tuple<BitsParam< 3>, uint32_t >,   //  12  1.5 uint64_t
		::std::tuple<BitsParam< 2>, uint64_t >,   //  16  (__m128i)
		::std::tuple<BitsParam< 3>, uint64_t >,   //  24  1.5 __m128i
		::std::tuple<BitsParam< 4>, uint64_t >,   //  32  (__m256i)
		::std::tuple<BitsParam< 5>, uint64_t >,   //  40  1.5 __m256i
		::std::tuple<BitsParam< 8>, uint64_t >    //  64  (512 bits)
> BitReverseTransformTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseTransformTest, BitReverseTransformTestTypes);

