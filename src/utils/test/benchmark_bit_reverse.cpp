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
#include <sstream>
#include <algorithm>

// include files to test
#include "utils/bitgroup_ops.hpp"

#include "utils/benchmark_utils.hpp"
#include "utils/test/bit_test_common.hpp"


//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bit reverse.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.

// NOTE if the gtest fixture class is missing the trailing semicolon, compile error about "expected initializer..."
template <unsigned char BITS_PER_GROUP>
class BitReverseBenchmarkHelper {
  public:

    uint8_t input[384] = {  0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                            0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                            0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10 };
    static constexpr size_t iters = 250000;
//    static constexpr size_t iters = 250;
};


//================== benchmark speed of SIMD_type + WORD_TPYE, fixing the bitgroup size and amount of data.

template <typename P>
class BitReverseWordBenchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<P::bitsPerGroup> helper;

    template <typename WORD_TYPE, uint8_t SIMD_TYPE, typename P2 = P, typename ::std::enable_if<(P2::bitsPerGroup < (sizeof(WORD_TYPE) * 8)), int>::type = 0>
    void word_test( std::string name ) {

	size_t count = 128;

        ::bliss::utils::bit_ops::bitgroup_ops<P2::bitsPerGroup, SIMD_TYPE> op;
	WORD_TYPE in[count];  // using "new" causes segv.  declare in thisway.
        for (size_t i = 0; i < count; ++i) {
          memcpy(&(in[i]), this->helper.input + i, sizeof(WORD_TYPE));
        }

        size_t max = BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters * count / sizeof(WORD_TYPE);

        //printf("max: %lu, count %lu, wordsize %lu\n", max, count / sizeof(WORD_TYPE), sizeof(WORD_TYPE));

        BL_TIMER_START(this->bitrev);
        // looping timer too slow.
//        BL_TIMER_LOOP_START(this->bitrev);
        for (size_t iter = 0; iter < max; ++iter) {
//          BL_TIMER_LOOP_RESUME(this->bitrev);
          in[(iter+67) % count] = op.reverse(in[iter % count]);  // + 67 so that we are on different cache lines. also, +64 means that after 128 iterations, we have the same locations - could be optimized out.
//          BL_TIMER_LOOP_PAUSE(this->bitrev);
        }  // else too large, so don't do the test.
//        BL_TIMER_LOOP_END(this->bitrev, name, BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters * 128);
        BL_TIMER_END(this->bitrev, name, BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters * 128);  // bytes


    }
    template <typename WORD_TYPE, uint8_t SIMD_TYPE, typename P2 = P, typename ::std::enable_if<(P2::bitsPerGroup >= (sizeof(WORD_TYPE) * 8)), int>::type = 0>
    void word_test( std::string name ) {
    }
  public:
#if BL_BENCHMARK == 1
    static BL_TIMER_INIT(bitrev);
#else
    BL_TIMER_INIT(bitrev); // does nothing
#endif
    static constexpr uint8_t bits = P::bitsPerGroup;

    static void TearDownTestCase() {
      BL_TIMER_REPORT(BitReverseWordBenchmark<P>::bitrev);
    }

};

#if BL_BENCHMARK == 1
template <typename P>
BL_TIMER_INIT(BitReverseWordBenchmark<P>::bitrev);
#endif
template <typename P>
constexpr uint8_t BitReverseWordBenchmark<P>::bits;


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseWordBenchmark);

TYPED_TEST_P(BitReverseWordBenchmark, reverse_word)
{

  this->template word_test<uint8_t , ::bliss::utils::bit_ops::BIT_REV_SEQ>("SEQ uint8");
  this->template word_test<uint16_t, ::bliss::utils::bit_ops::BIT_REV_SEQ>("SEQ uint16");
  this->template word_test<uint32_t, ::bliss::utils::bit_ops::BIT_REV_SEQ>("SEQ uint32");
  this->template word_test<uint64_t, ::bliss::utils::bit_ops::BIT_REV_SEQ>("SEQ uint64");
  this->template word_test<uint8_t , ::bliss::utils::bit_ops::BIT_REV_SWAR>("SWAR uint8");
  this->template word_test<uint16_t, ::bliss::utils::bit_ops::BIT_REV_SWAR>("SWAR uint16");
  this->template word_test<uint32_t, ::bliss::utils::bit_ops::BIT_REV_SWAR>("SWAR uint32");
  this->template word_test<uint64_t, ::bliss::utils::bit_ops::BIT_REV_SWAR>("SWAR uint64");

#ifdef __SSSE3__
  this->template word_test<__m128i, ::bliss::utils::bit_ops::BIT_REV_SSSE3>("SSSE m128i");
#endif

#ifdef __AVX2__
  this->template word_test<__m256i, ::bliss::utils::bit_ops::BIT_REV_AVX2>("AVX m256i");
#endif
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseWordBenchmark, reverse_word);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<
    ::bliss::utils::bit_ops::test::BitsParam< 1>,
    ::bliss::utils::bit_ops::test::BitsParam< 2>,
    ::bliss::utils::bit_ops::test::BitsParam< 3>,
    ::bliss::utils::bit_ops::test::BitsParam< 4>,
    ::bliss::utils::bit_ops::test::BitsParam< 8>,
    ::bliss::utils::bit_ops::test::BitsParam<16>,
    ::bliss::utils::bit_ops::test::BitsParam<32>,
    ::bliss::utils::bit_ops::test::BitsParam<64>,
    ::bliss::utils::bit_ops::test::BitsParam<128>
> BitReverseWordBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseWordBenchmark, BitReverseWordBenchmarkTypes);


//====================== fix bitgroup size and number of iters, check effect of array size (1 up to max) and simd type.
//   this is basically testing the remainder when converting a byte array.

template <typename P>
class BitReverseRemainderBenchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<P::bitsPerGroup> helper;


    template <typename WORD_TYPE, uint8_t SIMD_TYPE, typename P2 = P, typename ::std::enable_if<(P2::bitsPerGroup < (sizeof(WORD_TYPE) * 8)), int>::type = 0>
    void part_test( std::string name ) {

	::bliss::utils::bit_ops::bitgroup_ops<P2::bitsPerGroup, SIMD_TYPE> op;


  	size_t step = (SIMD_TYPE < 2) ? 8 : (SIMD_TYPE == 2) ? 16 : 32;
	size_t iters = BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters * 16;

	uint8_t BLISS_ALIGNED_ARRAY(out, 224, 32);
	memcpy(out, this->helper.input, 224);

	std::stringstream ss;	
	

	for (size_t i = 1; i <= step; ++i) {
	  if (((i * 8) % P2::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.
	  
	  ss.str(std::string());
  	  ss.clear();
	  ss << name << "_" << i;
	  	  
//          BL_TIMER_LOOP_START(this->bitrev);
          BL_TIMER_START(this->bitrev);


  	  for (size_t iter = 0; iter < iters; ++iter) {

//	    BL_TIMER_LOOP_RESUME(this->bitrev);
	    if ((P2::bitsPerGroup & (P2::bitsPerGroup - 1)) == 0)
              op.reverse(out + (iter + 67) % 134, out + iter % 134, i, 0);
	    else
              op.reverse(out + (iter + 67) % 134, out + iter % 134, i, iter % 8);

//            BL_TIMER_LOOP_PAUSE(this->bitrev);

          }
//          BL_TIMER_LOOP_END(this->bitrev, ss.str(), iters * i);
          BL_TIMER_END(this->bitrev, ss.str(), iters * i);
        }

    }
    template <typename WORD_TYPE, uint8_t SIMD_TYPE, typename P2 = P, typename ::std::enable_if<(P2::bitsPerGroup >= (sizeof(WORD_TYPE) * 8)), int>::type = 0>
    void part_test( std::string name ) {
    }

  public:
#if BL_BENCHMARK == 1
    static BL_TIMER_INIT(bitrev);
#else
    BL_TIMER_INIT(bitrev); // does nothing
#endif
    static constexpr uint8_t bits = P::bitsPerGroup;


    static void TearDownTestCase() {
      BL_TIMER_REPORT(BitReverseRemainderBenchmark<P>::bitrev);
    }

};

#if BL_BENCHMARK == 1
template <typename P>
BL_TIMER_INIT(BitReverseRemainderBenchmark<P>::bitrev);
#endif
template <typename P>
constexpr uint8_t BitReverseRemainderBenchmark<P>::bits;

// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseRemainderBenchmark);


TYPED_TEST_P(BitReverseRemainderBenchmark, reverse_remainder)
{
   this->template part_test<size_t, ::bliss::utils::bit_ops::BIT_REV_SEQ>("seq");
   this->template part_test<size_t, ::bliss::utils::bit_ops::BIT_REV_SWAR>("swar");
#ifdef __SSSE3__
   this->template part_test<__m128i, ::bliss::utils::bit_ops::BIT_REV_SSSE3>("ssse3");
#endif
#ifdef __AVX2__
   this->template part_test<__m256i, ::bliss::utils::bit_ops::BIT_REV_AVX2>("avx2");
#endif
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseRemainderBenchmark, reverse_remainder);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<
    ::bliss::utils::bit_ops::test::BitsParam< 1>,
    ::bliss::utils::bit_ops::test::BitsParam< 2>,
    ::bliss::utils::bit_ops::test::BitsParam< 3>,
    ::bliss::utils::bit_ops::test::BitsParam< 4>,
    ::bliss::utils::bit_ops::test::BitsParam< 8>,
    ::bliss::utils::bit_ops::test::BitsParam<16>,
    ::bliss::utils::bit_ops::test::BitsParam<32>,
    ::bliss::utils::bit_ops::test::BitsParam<64>,
    ::bliss::utils::bit_ops::test::BitsParam<128>
> BitReverseRemainderBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseRemainderBenchmark, BitReverseRemainderBenchmarkTypes);





//====================== fix bitgroup size and total number of bytes, check effect of maximum simd type.
//   this is testing converting the whole array.

template <typename P>
class BitReverseArrayBenchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<P::bitsPerGroup> helper;


    template <uint8_t MAX_SIMD_TYPE, typename P2 = P>
    void array_test( std::string name ) {
      uint8_t BLISS_ALIGNED_ARRAY(out, 384, 32);

      BL_TIMER_START(this->bitrev);
//      BL_TIMER_LOOP_START(this->bitrev);
      memcpy(out, this->helper.input, 384);

      for (size_t iter = 0; iter < BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters; ++iter) {

//        BL_TIMER_LOOP_RESUME(this->bitrev);
        bliss::utils::bit_ops::reverse<P2::bitsPerGroup, MAX_SIMD_TYPE>( out + (iter + 119) % 238, out + (iter % 238), 128);  // input from 0 to 256. output from 128-256, then 0 to 128
//        BL_TIMER_LOOP_PAUSE(this->bitrev);
      }
//      BL_TIMER_LOOP_END(this->bitrev, name, BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters * 128);
      BL_TIMER_END(this->bitrev, name, BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters * 128);
    }



  public:
#if BL_BENCHMARK == 1
    static BL_TIMER_INIT(bitrev);
#else
    BL_TIMER_INIT(bitrev); // does nothing
#endif
    static constexpr uint8_t bits = P::bitsPerGroup;

    static void TearDownTestCase() {
      BL_TIMER_REPORT(BitReverseArrayBenchmark<P>::bitrev);
    }

};

#if BL_BENCHMARK == 1
template <typename P>
BL_TIMER_INIT(BitReverseArrayBenchmark<P>::bitrev);
#endif
template <typename P>
constexpr uint8_t BitReverseArrayBenchmark<P>::bits;

// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseArrayBenchmark);


TYPED_TEST_P(BitReverseArrayBenchmark, reverse_short_array)
{
   this->template array_test<::bliss::utils::bit_ops::BIT_REV_SEQ>("seq");
   this->template array_test<::bliss::utils::bit_ops::BIT_REV_SWAR>("swar");
#ifdef __SSSE3__
   this->template array_test<::bliss::utils::bit_ops::BIT_REV_SSSE3>("ssse3");
#endif
#ifdef __AVX2__
   this->template array_test<::bliss::utils::bit_ops::BIT_REV_AVX2>("avx2");
#endif
}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseArrayBenchmark, reverse_short_array);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<
    ::bliss::utils::bit_ops::test::BitsParam< 1>,
    ::bliss::utils::bit_ops::test::BitsParam< 2>,
    ::bliss::utils::bit_ops::test::BitsParam< 3>,
    ::bliss::utils::bit_ops::test::BitsParam< 4>,
    ::bliss::utils::bit_ops::test::BitsParam< 8>,
    ::bliss::utils::bit_ops::test::BitsParam<16>,
    ::bliss::utils::bit_ops::test::BitsParam<32>
> BitReverseArrayBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseArrayBenchmark, BitReverseArrayBenchmarkTypes);




//====================== fix bitgroup size and total number of bytes, check effect of maximum simd type.
//   this is testing converting the whole array.

template <typename P>
class BitReverseFixedArrayBenchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<P::bitsPerGroup> helper;

    template <typename MAX_SIMD_TYPE, typename P2 = P>
    void array_test( std::string name ) {

    	constexpr size_t iters = BitReverseBenchmarkHelper<P2::bitsPerGroup>::iters;

      uint8_t data[10240][128];

      for (int i = 0; i < 10240; ++i) {
        memcpy(data[i], this->helper.input + (i % 256), 128);
      }

      constexpr uint16_t pad_bits = ((P2::bitsPerGroup & (P2::bitsPerGroup - 1)) == 0) ? 0 : ((128 * 8) % P2::bitsPerGroup);

      BL_TIMER_START(this->bitrev);
      for (size_t iter = 0; iter < iters; ++iter) {
        ::bliss::utils::bit_ops::template reverse<P2::bitsPerGroup, MAX_SIMD_TYPE, pad_bits>( data[iter % 10240], data[((iter + 1) % 10240)]);  // input from 0 to 256. output from 128-256, then 0 to 128
      }

      BL_TIMER_END(this->bitrev, name, iters * 128);
    }


  public:
#if BL_BENCHMARK == 1
    static BL_TIMER_INIT(bitrev);
#else
    BL_TIMER_INIT(bitrev); // does nothing
#endif
    static constexpr uint8_t bits = P::bitsPerGroup;

    static void TearDownTestCase() {
      BL_TIMER_REPORT(BitReverseFixedArrayBenchmark<P>::bitrev);
    }

};

#if BL_BENCHMARK == 1
template <typename P>
BL_TIMER_INIT(BitReverseFixedArrayBenchmark<P>::bitrev);
#endif
template <typename P>
constexpr uint8_t BitReverseFixedArrayBenchmark<P>::bits;

// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseFixedArrayBenchmark);


TYPED_TEST_P(BitReverseFixedArrayBenchmark, reverse_fixed_array)
{
   this->template array_test<::bliss::utils::bit_ops::BITREV_SEQ>("seq");

   this->template array_test<::bliss::utils::bit_ops::BITREV_SWAR>("swar");
#ifdef __SSSE3__
   this->template array_test<::bliss::utils::bit_ops::BITREV_SSSE3>("ssse3");
#endif
#ifdef __AVX2__
   this->template array_test<::bliss::utils::bit_ops::BITREV_AVX2>("avx2");
#endif


}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseFixedArrayBenchmark, reverse_fixed_array);


//////////////////// RUN the tests with different types.
typedef ::testing::Types<
    ::bliss::utils::bit_ops::test::BitsParam< 1>,
    ::bliss::utils::bit_ops::test::BitsParam< 2>,
    ::bliss::utils::bit_ops::test::BitsParam< 3>,
    ::bliss::utils::bit_ops::test::BitsParam< 4>,
    ::bliss::utils::bit_ops::test::BitsParam< 8>,
    ::bliss::utils::bit_ops::test::BitsParam<16>,
    ::bliss::utils::bit_ops::test::BitsParam<32>
> BitReverseFixedArrayBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseFixedArrayBenchmark, BitReverseFixedArrayBenchmarkTypes);








// special test for bit revrese within a byte.  compare to http://stackoverflow.com/questions/746171/best-algorithm-for-bit-reversal-from-msb-lsb-to-lsb-msb-in-c
// where best time is 0.05 for 400M bytes.  this code does 0.026s for 400M bytes (inrements of 1M) on a i5-4200U laptop.
// 0.39s in increments of 4M  (can't go higher without changing stack size (about 8M)
TEST(BitReverseInFixedByteArrayBenchmark, rev_byte)
{
	constexpr uint8_t bitsPerGroup = 1;

    BitReverseBenchmarkHelper<bitsPerGroup> helper;

    constexpr size_t iters = 400000000;
    constexpr size_t s = 4000000;

    uint8_t data[s];

    for (size_t i = 0; i < s; i+=100) {
      memcpy(data + i, helper.input + (i % 256), 100);
    }
//    size_t rem = iters % 100;
//    memcpy(data + (iters - iters % 128), helper.input, rem);

    BL_TIMER_INIT(bitrev); // does nothing

    BL_TIMER_START(bitrev);
    for (size_t iter = 0; iter < iters; iter += s) {
      ::bliss::utils::bit_ops::template reverse_bits_in_byte<1, ::bliss::utils::bit_ops::BITREV_SWAR, 0>(data, data);  // input from 0 to 256. output from 128-256, then 0 to 128
    }
    BL_TIMER_END(bitrev, "swar", iters);

    BL_TIMER_START(bitrev);
    for (size_t iter = 0; iter < iters; iter += s) {
      ::bliss::utils::bit_ops::template reverse_bits_in_byte<1, ::bliss::utils::bit_ops::BITREV_SSSE3, 0>(data, data);  // input from 0 to 256. output from 128-256, then 0 to 128
    }
    BL_TIMER_END(bitrev, "ssse3", iters);

    BL_TIMER_START(bitrev);
    for (size_t iter = 0; iter < iters; iter += s) {
      ::bliss::utils::bit_ops::template reverse_bits_in_byte<1, ::bliss::utils::bit_ops::BITREV_AVX2, 0>(data, data);  // input from 0 to 256. output from 128-256, then 0 to 128
    }
    BL_TIMER_END(bitrev, "avx2", iters);


    BL_TIMER_REPORT(bitrev);
}
