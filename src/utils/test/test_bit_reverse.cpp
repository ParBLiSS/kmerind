
// include google test
#include <gtest/gtest.h>

#include <random>
#include <array>
#include <cstdint>
#include <utility>
#include <iostream>

// include files to test
#include "utils/bitgroup_ops.hpp"

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



    template <unsigned char BITS = BITS_PER_GROUP, typename std::enable_if<(BITS >= 8), int>::type = 1>
    bool is_reverse(uint8_t const * orig, uint8_t const* rev, size_t len, uint8_t bit_offset = 0) {
      bool same = true;

      //std::cout << "u: " << std::hex << orig << " v: " << std::hex << rev << std::endl;


      // compare based on BITS size.
      // reverse bytes or larger, than just compare the bytes.
      const uint8_t* w = rev + len;
      const uint8_t* uu = orig;

      uint8_t bytesPerChar = BITS / 8;
      uint8_t max = len / bytesPerChar;

      bool local_same = false;
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

    template <unsigned char BITS = BITS_PER_GROUP, typename std::enable_if<(BITS < 8) && ((BITS & (BITS - 1)) != 0), int>::type = 1>
    bool is_reverse(uint8_t const * orig, uint8_t const* rev, size_t len, uint8_t bit_offset = 0) {
      bool same = true;

      if (len == 1) {  // single byte

        size_t rem = (8 - bit_offset) % BITS;
        unsigned char bits_mask = static_cast<unsigned char>(~(0xFF << BITS));
        //unsigned char rem_mask = static_cast<unsigned char>(~(0xFF >> rem));

        // check the last part is same.
        //same = (*orig & rem_mask) == (*rev & rem_mask);

        // check the reversed part
        // align test data so the reversed part is at the highest bits
        for (size_t i = bit_offset, j = (8 - bit_offset - BITS); i < (8 - rem); i += BITS, j -= BITS) {
          same &= ((*orig >> i) & bits_mask) == ((*rev >> j) & bits_mask);
        }

      } else {   // multibyte

        // use 2 bytes, overlapped by 1 byte in successive run, then all bits will be covered.

        uint8_t offset = bit_offset;  // offset from lowest.

        const uint8_t* w = rev + len - 1;  // doing it in pairs of 2.
        const uint8_t* uu = orig;

        uint16_t bits_mask = static_cast<uint16_t>(~(0xFFFF << BITS));

        bool local_same;

        for (uint8_t i = 0; i < len-1; ++i) {
          --w;

          size_t rem = (16 - offset) % BITS;  // remainder within the current short.

          for (uint8_t k = 16 - offset - BITS; offset < (16 - rem); offset += BITS, k -= BITS) {
            local_same = ((*(reinterpret_cast<const uint16_t*>(uu)) >> offset) & bits_mask) == ((*(reinterpret_cast<const uint16_t*>(w)) >> k) & bits_mask);

            if (!local_same) {
              std::cout << "diff @ " << std::dec << static_cast<size_t>(i) << ". off=" << static_cast<size_t>(offset) << ", k=" << static_cast<size_t>(k) << ":";
              std::cout << " full u=" << std::hex << *(reinterpret_cast<const uint16_t*>(uu)) << " full v=" << *(reinterpret_cast<const uint16_t*>(w));
              std::cout << " u=" << std::hex << static_cast<size_t>(((*(reinterpret_cast<const uint16_t*>(uu)) >> offset) & bits_mask));
              std::cout << " v=" << static_cast<size_t>(((*(reinterpret_cast<const uint16_t*>(w)) >> k) & bits_mask)) << std::endl;
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

      uint8_t orig_arr alignas(16) [16];
      uint8_t  rev_arr alignas(16) [16];

      _mm_store_si128((__m128i*)orig_arr, orig);
      _mm_store_si128((__m128i*)rev_arr, rev);


      return is_reverse(orig_arr,
                        rev_arr, 16, bit_offset);
    }
#endif

#ifdef __AVX2__

    bool is_reverse(__m256i const & orig, __m256i const & rev, uint8_t bit_offset = 0) {

      uint8_t orig_arr alignas(32) [32];
      uint8_t  rev_arr alignas(32) [32];

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

      EXPECT_TRUE(this->is_reverse(in, out));
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

      EXPECT_TRUE(this->is_reverse(in, out));
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

      EXPECT_TRUE(this->is_reverse(in, out));
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

      EXPECT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_short_array)
{

  TypeParam op;

  uint8_t out alignas(32) [32];

  unsigned int max = 8;

  for (unsigned int i = 1; i <= max; ++i ) {
    if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


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
          std::cout << "in: ";
          for (size_t j = 0; j < i; ++j) {
            std::cout << std::hex << static_cast<size_t>(this->helper.input[k + i - 1 - j]) << " ";

          }
          std::cout << std::endl;

          std::cout << "out: ";
          for (size_t j = 0; j < i; ++j) {
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
REGISTER_TYPED_TEST_CASE_P(BitReverseTest, reverse_uint8, reverse_uint16, reverse_uint32, reverse_uint64, reverse_short_array);


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

  uint8_t out alignas(32) [32];

  int max = 16;

  for (int i = 1; i <= max; ++i ) {
    if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


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

    EXPECT_TRUE(this->is_reverse(in, out));
  }  // else too large, so don't do the test.

}


TYPED_TEST_P(BitReverseAVX2Test, reverse_short_array)
{

  TypeParam op;

  uint8_t out alignas(32) [32];

  int max = 32;


  for (int i = 1; i <= max; ++i ) {
    if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


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

          EXPECT_TRUE(same);
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

        EXPECT_TRUE(same);
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
template <unsigned char Bits>
struct BitsParam { static constexpr unsigned char bitsPerGroup = Bits; };


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

  uint8_t out alignas(32) [128];

  unsigned int max = 128;

  for (unsigned int i = 1; i <= max; ++i ) {
    if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


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

