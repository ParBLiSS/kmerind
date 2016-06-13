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


#include "utils/test/bit_reverse_test_helper.hpp"

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of bit reverse.
//TESTS: for each, test different input (drawing from a 32 byte array),
//       different offsets, different bit group sizes, different word types, and different byte array lengths.
//TESTS: reverse entire array via multiplel SWAR, SSSE3, and AVX2 calls.





//
//
//template <typename T>
//class BitReverseFixedArrayTest : public ::testing::Test {
//  protected:
//
//    template <uint8_t BIT_GROUP_SIZE>
//    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len)  {
//      BitReverseTestHelper<BIT_GROUP_SIZE> helper;
//
//      return helper.is_reverse(out, in, len, 0);
////                               (len * 8) % BIT_GROUP_SIZE);  // if BITS = 3, offset needs to be calculated.
//    }
//
//    template <unsigned int BITS, uint8_t MAX_SIMD_TYPE,
//      typename data_type,
//      size_t data_size,
//      typename ::std::enable_if<((data_size * sizeof(data_type) * 8) <= BITS), int>::type = 1>
//    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {
//      if ((data_size * sizeof(data_type) * 8) == BITS) {
//        memcpy(out, in, data_size * sizeof(data_type)); return;
//      } else if ((data_size * sizeof(data_type) * 8) < BITS) {
//        return;
//      }
//    }
//    template <unsigned int BITS, uint8_t MAX_SIMD_TYPE,
//      typename data_type,
//      size_t data_size,
//      typename ::std::enable_if<((data_size * sizeof(data_type) * 8) > BITS), int>::type = 1>
//    void test(data_type (&in)[data_size], data_type (&out)[data_size]) {
//
//      BitReverseTestHelper<BITS> helper;
//      bool same = true;
//
//      for (unsigned int k = 0; k < (32 - data_size); ++k) {
////          printf("array size = %lu, sizeof(datatype) = %lu, bits = %u, SIMD = %u, k = %u\n", data_size, sizeof(data_type), BITS, MAX_SIMD_TYPE, k);
//
//        memcpy(in, helper.array + k, data_size * sizeof(data_type));
//
//        memset(out, 0, data_size * sizeof(data_type));
//
//        ::bliss::utils::bit_ops::template reverse<BITS, MAX_SIMD_TYPE>(out, in);
//
//        same =
//            this->is_reverse<BITS>(helper.array + k,
//                                   reinterpret_cast<uint8_t*>(out),
//                                   data_size * sizeof(data_type));
//
//        if (!same) {
//          std::cout << "in (MSB to LSB): ";
//          for (int64_t j = data_size - 1; j >= 0; --j) {
//            std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(in[j]) << " ";
//          }
//          std::cout << std::endl;
//
//          std::cout << "out (MSB to LSB): ";
//          for (int64_t j = data_size - 1; j >= 0; --j) {
//            std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(out[j]) << " ";
//
//          }
//          std::cout << std::endl;
//
//        }
//
//        ASSERT_TRUE(same);
//
//      }
//    }
//
//    template <uint8_t SIMDType, typename data_type, size_t data_size>
//    void run_tests() {
//    	  data_type in  [data_size];
//    	  data_type out [data_size];
//    	  // NEED TO SPECIFY data_type and data_size for icc.  not needed for clang or gcc.
//    	  switch (sizeof(data_type)) {
//    	    case 8:
//    	      test<32, SIMDType, data_type, data_size>(in, out);
//    	    case 4:
//    	      test<16, SIMDType, data_type, data_size>(in, out);
//    	    case 2:
//    	      test<8,  SIMDType, data_type, data_size>(in, out);
//    	    case 1:
//    	      test<4,  SIMDType, data_type, data_size>(in, out);
//    	      test<2,  SIMDType, data_type, data_size>(in, out);
//    	      test<1,  SIMDType, data_type, data_size>(in, out);
//    	      test<3,  SIMDType, data_type, data_size>(in, out);
//    	      break;
//    	    default:
//    	      break;
//    	  }
//
//    }
//};
//
//TYPED_TEST_CASE_P(BitReverseFixedArrayTest);
//
//TYPED_TEST_P(BitReverseFixedArrayTest, reverse_seq)
//{
//  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
//  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
//
//  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_SEQ, data_type, data_size>();
//}
//
//
//TYPED_TEST_P(BitReverseFixedArrayTest, reverse_swar)
//{
//  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
//  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
//
//  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_SWAR, data_type, data_size>();
//}
//
//#ifdef __SSSE3__
//TYPED_TEST_P(BitReverseFixedArrayTest, reverse_ssse3)
//{
//  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
//  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
//
//  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_SSSE3, data_type, data_size>();
//
//}
//#endif
//
//#ifdef __AVX2__
//TYPED_TEST_P(BitReverseFixedArrayTest, reverse_avx2)
//{
//  using data_type = typename ::std::tuple_element<1, TypeParam>::type;
//  constexpr size_t data_size = ::std::tuple_element<0, TypeParam>::type::bitsPerGroup;
//
//  this->template run_tests<::bliss::utils::bit_ops::BIT_REV_AVX2, data_type, data_size>();
//}
//#endif
//
//
//// now register the test cases
//REGISTER_TYPED_TEST_CASE_P(BitReverseFixedArrayTest,
//#ifdef __SSSE3__
//                           reverse_ssse3,
//#endif
//#ifdef __AVX2__
//                           reverse_avx2,
//#endif
//                           reverse_swar,
//                           reverse_seq);
//
//
////////////////////// RUN the tests with different types.
//
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
//> BitReverseFixedArrayTestTypes;
//INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseFixedArrayTest, BitReverseFixedArrayTestTypes);
//
//
//
//





template <typename T>
class BitReverseTransformTest : public ::testing::Test {
  protected:

    template <uint8_t BIT_GROUP_SIZE,  uint8_t pad_bits,
    	typename data_type,
    	size_t data_size >
    bool is_reverse(data_type (&in)[data_size], data_type (&out)[data_size])  {
      BitReverseTestHelper<BIT_GROUP_SIZE> helper;

      // shift output left by padbits
      using SIMD = bliss::utils::bit_ops::BITREV_AUTO_AGGRESSIVE<(sizeof(data_type) * data_size), bliss::utils::bit_ops::BITREV_SWAR>;

      data_type tmp[data_size];
      bliss::utils::bit_ops::left_shift<SIMD, pad_bits>(tmp, out);

      // clear the pad bits in input

      data_type tmp2[data_size];
      memcpy(tmp2, in, sizeof(data_type) * data_size);
      tmp2[data_size - 1] &= (::std::numeric_limits<data_type>::max() >> pad_bits);

//      std::cout << "tmp(in) (MSB to LSB): ";
//      for (int64_t j = data_size - 1; j >= 0; --j) {
//        std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(tmp2[j]) << " ";
//      }
//      std::cout << std::endl;
//
//      std::cout << "tmp(out) (MSB to LSB): ";
//      for (int64_t j = data_size - 1; j >= 0; --j) {
//        std::cout << std::hex << std::setw(2 * sizeof(data_type)) << std::setfill('0') << static_cast<size_t>(tmp[j]) << " ";
//
//      }
//      std::cout << std::endl;


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

          printf("array size = %lu, sizeof(datatype) = %lu, bits = %u, pad_bits = %u, SIMD = %u, input_offset = %u\n", data_size, sizeof(data_type), BITS, pad_bits, MAX_SIMD_TYPE::SIMDVal, k);
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
    	      test<32, SIMDType,  0, data_type, data_size>(in, out);
    	      test<32, SIMDType, 32, data_type, data_size>(in, out);
    	    case 4:
    	      test<16, SIMDType,  0, data_type, data_size>(in, out);
    	      test<16, SIMDType, 16, data_type, data_size>(in, out);
    	    case 2:
    	      test< 8, SIMDType,  0, data_type, data_size>(in, out);
    	      test< 8, SIMDType,  8, data_type, data_size>(in, out);
    	    case 1:
    	      test< 4, SIMDType,  0, data_type, data_size>(in, out);
    	      test< 4, SIMDType,  4, data_type, data_size>(in, out);
    	      test< 2, SIMDType,  0, data_type, data_size>(in, out);
    	      test< 2, SIMDType,  2, data_type, data_size>(in, out);
    	      test< 2, SIMDType,  6, data_type, data_size>(in, out);
    	      test< 1, SIMDType,  0, data_type, data_size>(in, out);
    	      test< 1, SIMDType,  1, data_type, data_size>(in, out);
    	      test< 1, SIMDType,  7, data_type, data_size>(in, out);
    	      test< 3, SIMDType,  ((sizeof(data_type) * data_size * 8) % 3 ), data_type, data_size>(in, out);
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




