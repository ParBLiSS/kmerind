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
#include "utils/test/bit_test_common.hpp"

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bitwise operations.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.



template <typename T>
class BitwiseOpsFixedArrayTest : public ::testing::Test {
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
    bool is_negated(data_type (&in)[data_size], data_type (&out)[data_size])  {

		bool result = true;

		for (size_t i = 0; i < data_size; ++i) {
			result &= (out[i] == static_cast<data_type>(~(in[i])));

			if (!result) {
				std::cout << "word " << std::dec << i << " out = " << std::hex <<
						static_cast<size_t>(out[i]) << " in = " <<
						static_cast<size_t>(in[i]) << " gold = " <<
						static_cast<size_t>(static_cast<data_type>(~(in[i]))) << std::endl;
			}
		}

		return result;
	}

	template <typename data_type, size_t data_size>
    bool is_anded(data_type (&lhs)[data_size], data_type (&rhs)[data_size], data_type (&out)[data_size])  {

		bool result = true;

		for (size_t i = 0; i < data_size; ++i) {
			result &= (out[i] == static_cast<data_type>(lhs[i] & rhs[i]));

			if (!result) {
				std::cout << "word " << std::dec << i << " out = " << std::hex <<
						static_cast<size_t>(out[i]) << " lhs = " <<
						static_cast<size_t>(lhs[i]) << " rhs = " <<
						static_cast<size_t>(rhs[i]) << " gold = " <<
						static_cast<size_t>(static_cast<data_type>(lhs[i] & rhs[i])) << std::endl;
			}

		}

		return result;
	}

	template <typename data_type, size_t data_size>
    bool is_ored(data_type (&lhs)[data_size], data_type (&rhs)[data_size], data_type (&out)[data_size])  {

		bool result = true;

		for (size_t i = 0; i < data_size; ++i) {
			result &= (out[i] == static_cast<data_type>(lhs[i] | rhs[i]));

			if (!result) {
				std::cout << "word " << std::dec << i << " out = " << std::hex <<
						static_cast<size_t>(out[i]) << " lhs = " <<
						static_cast<size_t>(lhs[i]) << " rhs = " <<
						static_cast<size_t>(rhs[i]) << " gold = " <<
						static_cast<size_t>(static_cast<data_type>(lhs[i] | rhs[i])) << std::endl;
			}
		}

		return result;
	}

	template <typename data_type, size_t data_size>
    bool is_xored(data_type (&lhs)[data_size], data_type (&rhs)[data_size], data_type (&out)[data_size])  {

		bool result = true;

		for (size_t i = 0; i < data_size; ++i) {
			result &= (out[i] == static_cast<data_type>(lhs[i] ^ rhs[i]));

			if (!result) {
				std::cout << "word " << std::dec << i << " out = " << std::hex <<
						static_cast<size_t>(out[i]) << " lhs = " <<
						static_cast<size_t>(lhs[i]) << " rhs = " <<
						static_cast<size_t>(rhs[i]) << " gold = " <<
						static_cast<size_t>(static_cast<data_type>(lhs[i] ^ rhs[i])) << std::endl;
			}

		}

		return result;
	}


    template <typename MAX_SIMD_TYPE,
      typename data_type,
      size_t data_size>
    void test(data_type (&lhs)[data_size], data_type (&rhs)[data_size], data_type (&out)[data_size], data_type (&out2)[data_size]) {

      static_assert(data_size < 32, "ERROR: testing only allows data size less than 32 elements");

      for (size_t k = 0; k < (32 - data_size - 1); ++k) {

        memcpy(lhs, this->array + k, data_size * sizeof(data_type));
        memcpy(rhs, this->array + k + 1, data_size * sizeof(data_type));
        memcpy(out2, this->array + k, data_size * sizeof(data_type));

        bool same = true;

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_not<MAX_SIMD_TYPE>(out, lhs);

        same = this->is_negated(lhs, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " negate in: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " negate out: ";
        	::bliss::utils::bit_ops::test::print(out);
        }
        ASSERT_TRUE(same);

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_and<MAX_SIMD_TYPE>(out, lhs, rhs);

        same = this->is_anded(lhs, rhs, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " and lhs: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " and rhs: ";
        	::bliss::utils::bit_ops::test::print(rhs);
        	std::cout << "iter " << std::dec << k << " and out: ";
        	::bliss::utils::bit_ops::test::print(out);
        }
        ASSERT_TRUE(same);



        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_or<MAX_SIMD_TYPE>(out, lhs, rhs);

        same = this->is_ored(lhs, rhs, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " or lhs: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " or rhs: ";
        	::bliss::utils::bit_ops::test::print(rhs);
        	std::cout << "iter " << std::dec << k << " or out: ";
        	::bliss::utils::bit_ops::test::print(out);
        }
        ASSERT_TRUE(same);


        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_xor<MAX_SIMD_TYPE>(out, lhs, rhs);

        same = this->is_xored(lhs, rhs, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " xor lhs: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " xor rhs: ";
        	::bliss::utils::bit_ops::test::print(rhs);
        	std::cout << "iter " << std::dec << k << " xor out: ";
        	::bliss::utils::bit_ops::test::print(out);
        }
        ASSERT_TRUE(same);



        memcpy(out2, this->array + k, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_not<MAX_SIMD_TYPE>(out2, out2);

        same = this->is_negated(lhs, out2);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " negate in: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " negate out: ";
        	::bliss::utils::bit_ops::test::print(out2);
        }
        ASSERT_TRUE(same);

        memcpy(out2, this->array + k, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_and<MAX_SIMD_TYPE>(out2, out2, rhs);

        same = this->is_anded(lhs, rhs, out2);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " and lhs: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " and rhs: ";
        	::bliss::utils::bit_ops::test::print(rhs);
        	std::cout << "iter " << std::dec << k << " and out: ";
        	::bliss::utils::bit_ops::test::print(out2);
        }
        ASSERT_TRUE(same);



        memcpy(out2, this->array + k, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_or<MAX_SIMD_TYPE>(out2, out2, rhs);

        same = this->is_ored(lhs, rhs, out2);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " or lhs: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " or rhs: ";
        	::bliss::utils::bit_ops::test::print(rhs);
        	std::cout << "iter " << std::dec << k << " or out: ";
        	::bliss::utils::bit_ops::test::print(out2);
        }
        ASSERT_TRUE(same);


        memcpy(out2, this->array + k, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template bit_xor<MAX_SIMD_TYPE>(out2, out2, rhs);

        same = this->is_xored(lhs, rhs, out2);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " xor lhs: ";
        	::bliss::utils::bit_ops::test::print(lhs);
        	std::cout << "iter " << std::dec << k << " xor rhs: ";
        	::bliss::utils::bit_ops::test::print(rhs);
        	std::cout << "iter " << std::dec << k << " xor out: ";
        	::bliss::utils::bit_ops::test::print(out2);
        }
        ASSERT_TRUE(same);

      }
    }

};

TYPED_TEST_CASE_P(BitwiseOpsFixedArrayTest);

TYPED_TEST_P(BitwiseOpsFixedArrayTest, bitwise_ops_seq)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
  data_type in  [data_size];
  data_type in2 [data_size];
  data_type out [data_size];
  data_type out2 [data_size];
  // NEED TO SPECIFY data_type and data_size for icc.  not needed for clang or gcc.
      this->template test<::bliss::utils::bit_ops::BITREV_SEQ, data_type, data_size>(in, in2, out, out2);
}


TYPED_TEST_P(BitwiseOpsFixedArrayTest, bitwise_ops_swar)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
  data_type in[data_size];
  data_type in2 [data_size];
  data_type out[data_size];
  data_type out2[data_size];

      this->template test<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size>(in, in2, out, out2);
}

#ifdef __SSSE3__
TYPED_TEST_P(BitwiseOpsFixedArrayTest, bitwise_ops_ssse3)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
  data_type in[data_size];
  data_type in2 [data_size];
  data_type out[data_size];
  data_type out2[data_size];

      this->template test<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size>(in, in2, out, out2);
}
#endif

#ifdef __AVX2__
TYPED_TEST_P(BitwiseOpsFixedArrayTest, bitwise_ops_avx2)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
  data_type in[data_size];
  data_type in2 [data_size];
  data_type out[data_size];
  data_type out2[data_size];

      this->template test<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size>(in, in2, out, out2);
}
#endif


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitwiseOpsFixedArrayTest,
#ifdef __SSSE3__
                           bitwise_ops_ssse3,
#endif
#ifdef __AVX2__
                           bitwise_ops_avx2,
#endif
                           bitwise_ops_swar,
                           bitwise_ops_seq);


//////////////////// RUN the tests with different types.

//typedef ::testing::Types<
//	    ::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint8_t >,    //  1  (uint8_t)
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint16_t >,   //  2  (uint16_t)
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 3>, uint8_t >,    //  3  1.5 uint16_t
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint32_t >,   //  4  (uint32_t)
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 3>, uint16_t >,   //  6  1.5 uint32_t
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint64_t >,   //  8   (uint64_t)
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 3>, uint32_t >,   //  12  1.5 uint64_t
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 2>, uint64_t >,   //  16  (__m128i)
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 3>, uint64_t >,   //  24  1.5 __m128i
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 4>, uint64_t >,   //  32  (__m256i)
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 5>, uint64_t >,   //  40  1.5 __m256i
//		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 8>, uint64_t >    //  64  (512 bits)
//> BitwiseOpsFixedArrayTestTypes;
//INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitwiseOpsFixedArrayTest, BitwiseOpsFixedArrayTestTypes);

