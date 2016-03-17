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
#include "utils/test/bitwise_ops_test_fixture.cpp"

typedef ::testing::Types<
		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 3>, uint16_t >,   //  6  1.5 uint32_t
		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 1>, uint64_t >,   //  8   (uint64_t)
		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 3>, uint32_t >,   //  12  1.5 uint64_t
		::std::tuple<::bliss::utils::bit_ops::test::BitsParam< 2>, uint64_t >   //  16  (__m128i)
> BitwiseOpsFixedArrayTestTypes2;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss_2, BitwiseOpsFixedArrayTest, BitwiseOpsFixedArrayTestTypes2);

