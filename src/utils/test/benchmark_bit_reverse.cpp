
// include google test
#include <gtest/gtest.h>
#include <utils/bit_reverse.hpp>

#include <random>
#include <array>
#include <cstdint>
#include <utility>
#include <iostream>

// include files to test
#include "utils/logging.h"


//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bit reverse.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.

// NOTE if the gtest fixture class is missing the trailing semicolon, compile error about "expected initializer..."
template <unsigned char BITS_PER_GROUP>
class BitReverseBenchmarkHelper {
  public:

    uint8_t input[32] = {  0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10 };
    size_t iters = 1000000;
};


template <typename T>
class BitReverseBenchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<T::bitsPerGroup> helper;

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseBenchmark);

TYPED_TEST_P(BitReverseBenchmark, reverse_uint8)
{
  using WORD_TYPE = uint8_t;

  TypeParam op;
  WORD_TYPE in[33 - sizeof(WORD_TYPE)];
  for (size_t i = 0; i <= (32 - sizeof(WORD_TYPE)); ++i) {
	  memcpy(&(in[i]), this->helper.input + i, sizeof(WORD_TYPE));
  }

    if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
      for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {


        for (size_t iter = 0; iter < this->helper.iters; ++iter) {
          in[(k+1)%(33-sizeof(WORD_TYPE))] = op.operator()(in[k]);

      }
    }  // else too large, so don't do the test.
  }
	std::cout << static_cast<size_t>(in[0]) << std::endl;
}

TYPED_TEST_P(BitReverseBenchmark, reverse_uint16)
{
  using WORD_TYPE = uint16_t;

  TypeParam op;
  WORD_TYPE in[33 - sizeof(WORD_TYPE)];
  for (size_t i = 0; i <= (32 - sizeof(WORD_TYPE)); ++i) {
	  memcpy(&(in[i]), this->helper.input + i, sizeof(WORD_TYPE));
  }

    if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
      for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {


        for (size_t iter = 0; iter < this->helper.iters; ++iter) {
          in[(k+1)%(33-sizeof(WORD_TYPE))] = op.operator()(in[k]);

      }
    }  // else too large, so don't do the test.
  }
	std::cout << static_cast<size_t>(in[0]) << std::endl;
}

TYPED_TEST_P(BitReverseBenchmark, reverse_uint32)
{
  using WORD_TYPE = uint32_t;

  TypeParam op;
  WORD_TYPE in[33 - sizeof(WORD_TYPE)];
  for (size_t i = 0; i <= (32 - sizeof(WORD_TYPE)); ++i) {
	  memcpy(&(in[i]), this->helper.input + i, sizeof(WORD_TYPE));
  }

    if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
      for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {


        for (size_t iter = 0; iter < this->helper.iters; ++iter) {
          in[(k+1)%(33-sizeof(WORD_TYPE))] = op.operator()(in[k]);

      }
    }  // else too large, so don't do the test.
  }
	std::cout << static_cast<size_t>(in[0]) << std::endl;
}


TYPED_TEST_P(BitReverseBenchmark, reverse_uint64)
{
  using WORD_TYPE = uint64_t;

  TypeParam op;
  WORD_TYPE in[33 - sizeof(WORD_TYPE)];
  for (size_t i = 0; i <= (32 - sizeof(WORD_TYPE)); ++i) {
	  memcpy(&(in[i]), this->helper.input + i, sizeof(WORD_TYPE));
  }

    if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
      for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {


        for (size_t iter = 0; iter < this->helper.iters; ++iter) {
          in[(k+1)%(33-sizeof(WORD_TYPE))] = op.operator()(in[k]);

      }
    }  // else too large, so don't do the test.
  }
	std::cout << static_cast<size_t>(in[0]) << std::endl;
}

TYPED_TEST_P(BitReverseBenchmark, reverse_short_array)
{
  TypeParam op;

  uint8_t out[32] alignas(32);
  memcpy(out, this->helper.input, 32);

  int max = 8;

  if (TypeParam::bitsPerGroup == 3) {
	  size_t iter = 0;

	  while (iter < this->helper.iters) {
		  for (int i = 1; i <= max; ++i ) {
			  if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.

			  for (int k = 0; k <= (32 - i); ++k) {
				  for (int j = 0; j < 8; ++j) {

					  op.operator()(out, out + k, i, j);
					  ++iter;
//					  if (iter % 1000 == 0) printf(".");
					  if (iter >= this->helper.iters) break;
				  }
				  if (iter >= this->helper.iters) break;

			  }
			  if (iter >= this->helper.iters) break;

		  }
	  }
	  for (int i = 0; i < 32; ++i) {
		  std::cout << static_cast<size_t>(out[31 - i]) << " ";
	  }
	  std::cout << std::endl;
  } else {
	  size_t iter = 0;

	  while (iter < this->helper.iters) {
		  for (int i = 1; i <= max; ++i ) {
			  if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.

			  for (int k = 0; k <= (32 - i); ++k) {

				  op.operator()(out, out + k, i, 0);
				  ++iter;
//				  if (iter % 1000 == 0) printf(".");
				  if (iter >= this->helper.iters) break;
			  }
			  if (iter >= this->helper.iters) break;
		  }
	  }
	  for (int i = 0; i < 32; ++i) {
		  std::cout << static_cast<size_t>(out[31 - i]) << " ";
	  }
	  std::cout << std::endl;
  }

}



// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseBenchmark, reverse_uint8, reverse_uint16, reverse_uint32, reverse_uint64, reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bit_reverse< 1, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
     ::bliss::utils::bit_ops::bit_reverse< 2, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
      ::bliss::utils::bit_ops::bit_reverse< 3, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
       ::bliss::utils::bit_ops::bit_reverse< 4, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
        ::bliss::utils::bit_ops::bit_reverse< 8, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
         ::bliss::utils::bit_ops::bit_reverse<16, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
          ::bliss::utils::bit_ops::bit_reverse<32, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
           ::bliss::utils::bit_ops::bit_reverse< 1, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
            ::bliss::utils::bit_ops::bit_reverse< 2, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
             ::bliss::utils::bit_ops::bit_reverse< 3, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
              ::bliss::utils::bit_ops::bit_reverse< 4, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
               ::bliss::utils::bit_ops::bit_reverse< 8, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
                ::bliss::utils::bit_ops::bit_reverse<16, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
                 ::bliss::utils::bit_ops::bit_reverse<32, ::bliss::utils::bit_ops::BIT_REV_SWAR>
> BitReverseBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseBenchmark, BitReverseBenchmarkTypes);



#ifdef __SSSE3__

template <typename T>
class BitReverseSSSEBenchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<T::bitsPerGroup> helper;


};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseSSSEBenchmark);

TYPED_TEST_P(BitReverseSSSEBenchmark, reverse_m128i)
{
  TypeParam op;

  __m128i in[17];

  for (int i = 0; i <= 16; ++i) {
	  in[i] = _mm_loadu_si128((__m128i*)(this->helper.input + i));

  }

    if (TypeParam::bitsPerGroup < 128 ) {
      for (size_t k = 0; k <= 16; ++k) {


        for (size_t iter = 0; iter < this->helper.iters; ++iter) {
          in[(k+1)%17] = op.operator()(in[k]);

      }
    }  // else too large, so don't do the test.
  }
	std::cout << op.toString(in[0]) << std::endl;
}


TYPED_TEST_P(BitReverseSSSEBenchmark, reverse_short_array)
{

  TypeParam op;

  uint8_t out[32] alignas(32);
  memcpy(out, this->helper.input, 32);

  int max = 16;

  if (TypeParam::bitsPerGroup == 3) {
	  size_t iter = 0;

	  while (iter < this->helper.iters) {

		  for (int i = 1; i <= max; ++i ) {
			  if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.

			  for (int k = 0; k <= (32 - i); ++k) {
				  for (int j = 0; j < 8; ++j) {

					  op.operator()(out, out + k, i, j);
					  ++iter;
//					  if (iter % 1000 == 0) printf(".");
					  if (iter >= this->helper.iters) break;
				  }
				  if (iter >= this->helper.iters) break;

			  }
			  if (iter >= this->helper.iters) break;

		  }
	  }
	  for (int i = 0; i < 32; ++i) {
		  std::cout << static_cast<size_t>(out[31 - i]) << " ";
	  }
	  std::cout << std::endl;
  } else {
	  size_t iter = 0;

	  while (iter < this->helper.iters) {

		  for (int i = 1; i <= max; ++i ) {
			  if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.

			  for (int k = 0; k <= (32 - i); ++k) {

				  op.operator()(out, out + k, i, 0);
				  ++iter;
//				  if (iter % 1000 == 0) printf(".");
				  if (iter >= this->helper.iters) break;
			  }
			  if (iter >= this->helper.iters) break;
		  }
	  }
	  for (int i = 0; i < 32; ++i) {
		  std::cout << static_cast<size_t>(out[31 - i]) << " ";
	  }
	  std::cout << std::endl;
  }

}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseSSSEBenchmark, reverse_m128i, reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bit_reverse< 1, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
     ::bliss::utils::bit_ops::bit_reverse< 2, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
      ::bliss::utils::bit_ops::bit_reverse< 3, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
       ::bliss::utils::bit_ops::bit_reverse< 4, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
        ::bliss::utils::bit_ops::bit_reverse< 8, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
         ::bliss::utils::bit_ops::bit_reverse<16, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
          ::bliss::utils::bit_ops::bit_reverse<32, ::bliss::utils::bit_ops::BIT_REV_SSSE3> ,
           ::bliss::utils::bit_ops::bit_reverse<64, ::bliss::utils::bit_ops::BIT_REV_SSSE3>
> BitReverseSSSEBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseSSSEBenchmark, BitReverseSSSEBenchmarkTypes);

#endif



#ifdef __AVX2__
template <typename T>
class BitReverseAVX2Benchmark : public ::testing::Test {
  protected:

    BitReverseBenchmarkHelper<T::bitsPerGroup> helper;

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseAVX2Benchmark);

TYPED_TEST_P(BitReverseAVX2Benchmark, reverse_m256i)
{
	TypeParam op;
  __m256i in = _mm256_loadu_si256((__m256i*)(this->helper.input));

    if (TypeParam::bitsPerGroup < 256 ) {

      for (size_t iter = 0; iter < this->helper.iters; ++iter) {
        in = op.operator()(in);

    }  // else too large, so don't do the test.
  }
	std::cout << op.toString(in) << std::endl;
}


TYPED_TEST_P(BitReverseAVX2Benchmark, reverse_short_array)
{

  TypeParam op;

  uint8_t out[32] alignas(32);
  memcpy(out, this->helper.input, 32);

  int max = 32;

  if (TypeParam::bitsPerGroup == 3) {
	  size_t iter = 0;

	  while (iter < this->helper.iters) {

		  for (int i = 1; i <= max; ++i ) {
			  if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.

			  for (int k = 0; k <= (32 - i); ++k) {
				  for (int j = 0; j < 8; ++j) {

					  op.operator()(out, out + k, i, j);
					  ++iter;

//					  if (iter % 1000 == 0) printf(".");
					  if (iter >= this->helper.iters) break;
				  }
				  if (iter >= this->helper.iters) break;

			  }
			  if (iter >= this->helper.iters) break;

		  }
	  }
	  for (int i = 0; i < 32; ++i) {
		  std::cout << static_cast<size_t>(out[31 - i]) << " ";
	  }
	  std::cout << std::endl;
  } else {
	  size_t iter = 0;

	  while (iter < this->helper.iters) {

		  for (int i = 1; i <= max; ++i ) {
			  if ((i % ((TypeParam::bitsPerGroup + 7) / 8)) > 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.

			  for (int k = 0; k <= (32 - i); ++k) {

				  op.operator()(out, out + k, i, 0);
				  ++iter;
//				  if (iter % 1000 == 0) printf(".");
				  if (iter >= this->helper.iters) break;
			  }
			  if (iter >= this->helper.iters) break;
		  }
	  }
	  for (int i = 0; i < 32; ++i) {
		  std::cout << static_cast<size_t>(out[31 - i]) << " ";
	  }
	  std::cout << std::endl;
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseAVX2Benchmark, reverse_m256i, reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bit_reverse< 1, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
     ::bliss::utils::bit_ops::bit_reverse< 2, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
      ::bliss::utils::bit_ops::bit_reverse< 3, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
       ::bliss::utils::bit_ops::bit_reverse< 4, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
        ::bliss::utils::bit_ops::bit_reverse< 8, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
         ::bliss::utils::bit_ops::bit_reverse<16, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
          ::bliss::utils::bit_ops::bit_reverse<32, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
           ::bliss::utils::bit_ops::bit_reverse<64, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
            ::bliss::utils::bit_ops::bit_reverse<128, ::bliss::utils::bit_ops::BIT_REV_AVX2>
> BitReverseAVX2BenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseAVX2Benchmark, BitReverseAVX2BenchmarkTypes);

#endif





