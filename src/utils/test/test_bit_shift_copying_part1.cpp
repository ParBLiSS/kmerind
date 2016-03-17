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

#include "utils/test/bit_shift_test_fixture_copying.cpp"


//////////////////// RUN the tests with different types.

typedef ::testing::Types<   // select power of 2 up to 512 bits, fill in-between so that each 1/2 machine word type is represented
	    ::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint8_t >,    //  1  (uint8_t)
		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint16_t >   //  2  (uint16_t)
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
> BitShiftFixedArrayTestTypes1;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss_1, BitShiftFixedArrayCopyingTest, BitShiftFixedArrayTestTypes1);

