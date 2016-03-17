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



template <typename T>
class BitShiftFixedArrayCopyingTest : public ::testing::Test {
  protected:


    uint8_t array[128] = { 0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10,
                           0x1F,0x3E,0x5D,0x7C,0x9B,0xBA,0xD9,0xF8,0xE7,0xC6,0xA5,0x84,0x63,0x42,0x21,0x00,
                           0x0F,0x2E,0x4D,0x6C,0x8B,0xAA,0xC9,0xE8,0xF7,0xD6,0xB5,0x94,0x73,0x52,0x31,0x10};



    template <unsigned int BITS, typename MAX_SIMD_TYPE,
      typename data_type,
      size_t data_size>
    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {

      static_assert(data_size < 32, "ERROR: testing only allows data size less than 32 elements");

      for (size_t k = 0; k < (32 - data_size); ++k) {

        memcpy(in, this->array + k, data_size * sizeof(data_type));
        bool same = true;

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template right_shift<MAX_SIMD_TYPE, BITS>(out, in);

        same = ::bliss::utils::bit_ops::test::is_right_shifted<BITS>(in, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " right " << std::dec << BITS << " in: ";
        	::bliss::utils::bit_ops::test::print(in);
        	std::cout << "iter " << std::dec << k << " right " << std::dec << BITS << " out: ";
        	::bliss::utils::bit_ops::test::print(out);
        }
        ASSERT_TRUE(same);

        memset(out, 0, data_size * sizeof(data_type));

        ::bliss::utils::bit_ops::template left_shift<MAX_SIMD_TYPE, BITS>(out, in);

        same = ::bliss::utils::bit_ops::test::is_left_shifted<BITS>(in, out);

        if (!same) {
        	std::cout << "SIMD type " << std::dec << MAX_SIMD_TYPE::SIMDVal << std::endl;
        	std::cout << "iter " << std::dec << k << " left " << std::dec << BITS << " in: ";
        	::bliss::utils::bit_ops::test::print(in);
        	std::cout << "iter " << std::dec << k << " left " << std::dec << BITS << " out: ";
        	::bliss::utils::bit_ops::test::print(out);
        }
        ASSERT_TRUE(same);

      }
    }

    template <typename SIMDType, typename data_type, size_t data_size,
     typename ::std::enable_if<(sizeof(data_type) >= 8), int>::type = 1>
    void run_tests_8() {

    	  data_type in[data_size];
    	  data_type out[data_size];

		  test<63, SIMDType, data_type, data_size>(in, out);
		  test<48, SIMDType, data_type, data_size>(in, out);
		  test<33, SIMDType, data_type, data_size>(in, out);
		  test<32, SIMDType, data_type, data_size>(in, out);
    }
    template <typename SIMDType, typename data_type, size_t data_size,
     typename ::std::enable_if<(sizeof(data_type) < 8), int>::type = 1>
    void run_tests_8() { }


    template <typename SIMDType, typename data_type, size_t data_size,
    	typename ::std::enable_if<sizeof(data_type) >= 4, int>::type = 1>
    void run_tests_4() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	  	  test<31, SIMDType, data_type, data_size>(in, out);
    	      test<24, SIMDType, data_type, data_size>(in, out);
    	      test<17, SIMDType, data_type, data_size>(in, out);
    	      test<16, SIMDType, data_type, data_size>(in, out);
    }
    template <typename SIMDType, typename data_type, size_t data_size,
    	typename ::std::enable_if<sizeof(data_type) < 4, int>::type = 1>
    void run_tests_4() { }

    template <typename SIMDType, typename data_type, size_t data_size,
	typename ::std::enable_if<sizeof(data_type) >= 2, int>::type = 1>
    void run_tests_2() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	      test<15, SIMDType, data_type, data_size>(in, out);
    	      test<12, SIMDType, data_type, data_size>(in, out);
    	      test< 9, SIMDType, data_type, data_size>(in, out);
    	      test< 8, SIMDType, data_type, data_size>(in, out);
    }
    template <typename SIMDType, typename data_type, size_t data_size,
	typename ::std::enable_if<sizeof(data_type) < 2, int>::type = 1>
    void run_tests_2() {}

    template <typename SIMDType, typename data_type, size_t data_size>
    void run_tests_1() {

    	  data_type in[data_size];
    	  data_type out[data_size];

    	      test< 7, SIMDType, data_type, data_size>(in, out);
    	      test< 6, SIMDType, data_type, data_size>(in, out);
    	      test< 5, SIMDType, data_type, data_size>(in, out);
    	      test< 4, SIMDType, data_type, data_size>(in, out);
    	      test< 3, SIMDType, data_type, data_size>(in, out);
    	      test< 2, SIMDType, data_type, data_size>(in, out);
    	      test< 1, SIMDType, data_type, data_size>(in, out);
    	      test< 0, SIMDType, data_type, data_size>(in, out);
    }
};

TYPED_TEST_CASE_P(BitShiftFixedArrayCopyingTest);


TYPED_TEST_P(BitShiftFixedArrayCopyingTest, shift_swar)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests_1<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size>();
  this->template run_tests_2<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size>();
  this->template run_tests_4<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size>();
  this->template run_tests_8<::bliss::utils::bit_ops::BITREV_SWAR, data_type, data_size>();
}

#ifdef __SSSE3__
TYPED_TEST_P(BitShiftFixedArrayCopyingTest, shift_ssse3)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests_1<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size>();
  this->template run_tests_2<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size>();
  this->template run_tests_4<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size>();
  this->template run_tests_8<::bliss::utils::bit_ops::BITREV_SSSE3, data_type, data_size>();
}
#endif

#ifdef __AVX2__
TYPED_TEST_P(BitShiftFixedArrayCopyingTest, shift_avx2)
{
  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;

  this->template run_tests_1<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size>();
  this->template run_tests_2<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size>();
  this->template run_tests_4<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size>();
  this->template run_tests_8<::bliss::utils::bit_ops::BITREV_AVX2, data_type, data_size>();
}
#endif


// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitShiftFixedArrayCopyingTest,
#ifdef __SSSE3__
                           shift_ssse3,
#endif
#ifdef __AVX2__
                           shift_avx2,
#endif
                           shift_swar
                           );

