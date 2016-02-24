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

// include files to test
#include "utils/bitgroup_ops.hpp"

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bit reverse.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.

template <unsigned char Bits>
struct BitsParam { static constexpr unsigned char bitsPerGroup = Bits; };







template <typename T>
class BitShiftFixedArrayTest : public ::testing::Test {
  protected:


    uint8_t array[128] = { 0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10};

	template <typename data_type, size_t data_size>
    void print(data_type (&w)[data_size]) {
      std::cout << "data type size " << std::dec << sizeof(data_type) << " len " << data_size << ": ";
      for (int k = data_size -1 ; k >= 0; --k) {
        std::cout << std::hex << static_cast<size_t>(w[k]) << " ";
      }
      std::cout << std::endl;
    }


	template <uint16_t shift, typename data_type, size_t data_size>
    bool is_right_shifted(data_type (&in)[data_size], data_type (&out)[data_size])  {

		if (shift == 0) {
			return std::equal(in, in + data_size, out);
		}

		data_type prev = 0;
		uint16_t inv_shift = sizeof(data_type) * 8 - shift;

		bool result = true;

		for (int i = data_size - 1; i >= 0; --i) {
			result &= (out[i] == static_cast<data_type>((in[i] >> shift) | prev));

			if (!result) {
				std::cout << "word " << std::dec << i << " out = " << std::hex <<
						static_cast<size_t>(out[i]) << " in = " <<
						static_cast<size_t>(in[i]) << " gold = " <<
						static_cast<size_t>((in[i] >> shift) | prev) << std::endl;
			}


			prev = in[i] << inv_shift;
		}

		return result;
	}

	template <uint16_t shift, typename data_type, size_t data_size>
    bool is_left_shifted(data_type (&in)[data_size], data_type (&out)[data_size])  {

		if (shift == 0) {
			return std::equal(in, in + data_size, out);
		}

		data_type prev = 0;
		uint16_t inv_shift = sizeof(data_type) * 8 - shift;

		bool result = true;

		for (size_t i = 0; i < data_size; ++i) {
			result &= (out[i] == static_cast<data_type>((in[i] << shift) | prev));

			if (!result) {
				std::cout << "word " << std::dec << i << " out = " << std::hex <<
						static_cast<size_t>(out[i]) << " in = " <<
						static_cast<size_t>(in[i]) << " gold = " <<
						static_cast<size_t>((in[i] << shift) | prev) << std::endl;
			}

			prev = in[i] >> inv_shift;
		}

		return result;
	}


    template <unsigned int BITS, typename MAX_SIMD_TYPE, bool inplace,
      typename data_type,
      size_t data_size, typename std::enable_if<inplace, int>::type = 1>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {

      for (unsigned int k = 0; k < (32 - data_size); ++k) {

        memcpy(in, this->array + k, data_size * sizeof(data_type));
        bool same = true;


        memcpy(out, in, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template right_shift<MAX_SIMD_TYPE, BITS>(out, out);

        same = this->is_right_shifted<BITS>(in, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " right inplace " << std::dec << BITS << " in: ";
        	print(in);
        	std::cout << "iter " << std::dec << k << " right inplace " << std::dec << BITS << " out: ";
        	print(out);
        }
        ASSERT_TRUE(same);

        memcpy(out, in, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template left_shift<MAX_SIMD_TYPE, BITS>(out, out);

        same = this->is_left_shifted<BITS>(in, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " left inplace " << std::dec << BITS << " in: ";
        	print(in);
        	std::cout << "iter " << std::dec << k << " left inplace " << std::dec << BITS << " out: ";
        	print(out);
        }

        ASSERT_TRUE(same);

      }
    }

    template <unsigned int BITS, typename MAX_SIMD_TYPE, bool inplace,
      typename data_type,
      size_t data_size, typename std::enable_if<!inplace, int>::type = 1>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {

      for (unsigned int k = 0; k < (32 - data_size); ++k) {

        memcpy(in, this->array + k, data_size * sizeof(data_type));
        bool same = true;

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template right_shift<MAX_SIMD_TYPE, BITS>(out, in);

        same = this->is_right_shifted<BITS>(in, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " right " << std::dec << BITS << " in: ";
        	print(in);
        	std::cout << "iter " << std::dec << k << " right " << std::dec << BITS << " out: ";
        	print(out);
        }
        ASSERT_TRUE(same);

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template left_shift<MAX_SIMD_TYPE, BITS>(out, in);

        same = this->is_left_shifted<BITS>(in, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " left " << std::dec << BITS << " in: ";
        	print(in);
        	std::cout << "iter " << std::dec << k << " left " << std::dec << BITS << " out: ";
        	print(out);
        }
        ASSERT_TRUE(same);

      }
    }

    template <typename SIMDType, typename data_type, size_t data_size, bool inplace,
     typename ::std::enable_if<(sizeof(data_type) >= 8), int>::type = 1>
    void run_tests() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	      test<63, SIMDType, inplace>(in, out);
    	      test<48, SIMDType, inplace>(in, out);
    	      test<33, SIMDType, inplace>(in, out);
    	      test<32, SIMDType, inplace>(in, out);
    	  	  test<31, SIMDType, inplace>(in, out);
    	      test<24, SIMDType, inplace>(in, out);
    	      test<17, SIMDType, inplace>(in, out);
    	      test<16, SIMDType, inplace>(in, out);
    	      test<15, SIMDType, inplace>(in, out);
    	      test<12, SIMDType, inplace>(in, out);
    	      test< 9, SIMDType, inplace>(in, out);
    	      test< 8, SIMDType, inplace>(in, out);
    	      test< 7, SIMDType, inplace>(in, out);
    	      test< 6, SIMDType, inplace>(in, out);
    	      test< 5, SIMDType, inplace>(in, out);
    	      test< 4, SIMDType, inplace>(in, out);
    	      test< 3, SIMDType, inplace>(in, out);
    	      test< 2, SIMDType, inplace>(in, out);
    	      test< 1, SIMDType, inplace>(in, out);
    	      test< 0, SIMDType, inplace>(in, out);
    }
    template <typename SIMDType, typename data_type, size_t data_size, bool inplace,
    	typename ::std::enable_if<(sizeof(data_type) >= 4) && (sizeof(data_type) < 8), int>::type = 1>
    void run_tests() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	  	  test<31, SIMDType, inplace>(in, out);
    	      test<24, SIMDType, inplace>(in, out);
    	      test<17, SIMDType, inplace>(in, out);
    	      test<16, SIMDType, inplace>(in, out);
    	      test<15, SIMDType, inplace>(in, out);
    	      test<12, SIMDType, inplace>(in, out);
    	      test< 9, SIMDType, inplace>(in, out);
    	      test< 8, SIMDType, inplace>(in, out);
    	      test< 7, SIMDType, inplace>(in, out);
    	      test< 6, SIMDType, inplace>(in, out);
    	      test< 5, SIMDType, inplace>(in, out);
    	      test< 4, SIMDType, inplace>(in, out);
    	      test< 3, SIMDType, inplace>(in, out);
    	      test< 2, SIMDType, inplace>(in, out);
    	      test< 1, SIMDType, inplace>(in, out);
    	      test< 0, SIMDType, inplace>(in, out);
    }
    template <typename SIMDType, typename data_type, size_t data_size, bool inplace,
	typename ::std::enable_if<(sizeof(data_type) >= 2) && (sizeof(data_type) < 4), int>::type = 1>
    void run_tests() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	      test<15, SIMDType, inplace>(in, out);
    	      test<12, SIMDType, inplace>(in, out);
    	      test< 9, SIMDType, inplace>(in, out);
    	      test< 8, SIMDType, inplace>(in, out);
    	      test< 7, SIMDType, inplace>(in, out);
    	      test< 6, SIMDType, inplace>(in, out);
    	      test< 5, SIMDType, inplace>(in, out);
    	      test< 4, SIMDType, inplace>(in, out);
    	      test< 3, SIMDType, inplace>(in, out);
    	      test< 2, SIMDType, inplace>(in, out);
    	      test< 1, SIMDType, inplace>(in, out);
    	      test< 0, SIMDType, inplace>(in, out);
    }
    template <typename SIMDType, typename data_type, size_t data_size, bool inplace,
	typename ::std::enable_if<(sizeof(data_type) < 2), int>::type = 1>
    void run_tests() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	      test< 7, SIMDType, inplace>(in, out);
    	      test< 6, SIMDType, inplace>(in, out);
    	      test< 5, SIMDType, inplace>(in, out);
    	      test< 4, SIMDType, inplace>(in, out);
    	      test< 3, SIMDType, inplace>(in, out);
    	      test< 2, SIMDType, inplace>(in, out);
    	      test< 1, SIMDType, inplace>(in, out);
    	      test< 0, SIMDType, inplace>(in, out);
    }
};

TYPED_TEST_CASE_P(BitShiftFixedArrayTest);


TYPED_TEST_P(BitShiftFixedArrayTest, shift_swar)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size, false>();
}

#ifdef __SSSE3__
TYPED_TEST_P(BitShiftFixedArrayTest, shift_ssse3)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size, false>();
}
#endif

#ifdef __AVX2__
TYPED_TEST_P(BitShiftFixedArrayTest, shift_avx2)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size, false>();
}
#endif


TYPED_TEST_P(BitShiftFixedArrayTest, shift_swar_inplace)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size, true>();
}

#ifdef __SSSE3__
TYPED_TEST_P(BitShiftFixedArrayTest, shift_ssse3_inplace)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size, true>();
}
#endif

#ifdef __AVX2__
TYPED_TEST_P(BitShiftFixedArrayTest, shift_avx2_inplace)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size, true>();
}
#endif


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitShiftFixedArrayTest,
#ifdef __SSSE3__
                           shift_ssse3,
                           shift_ssse3_inplace,
#endif
#ifdef __AVX2__
                           shift_avx2,
                           shift_avx2_inplace,
#endif
                           shift_swar,
                           shift_swar_inplace);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<   // select power of 2 up to 512 bits, fill in-between so that each 1/2 machine word type is represented
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
> BitShiftFixedArrayTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitShiftFixedArrayTest, BitShiftFixedArrayTestTypes);

