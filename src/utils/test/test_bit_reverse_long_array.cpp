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
    ::bliss::utils::bit_ops::test::BitsParam< 1> ,
    ::bliss::utils::bit_ops::test::BitsParam< 2> ,
    ::bliss::utils::bit_ops::test::BitsParam< 3> ,
    ::bliss::utils::bit_ops::test::BitsParam< 4> ,
    ::bliss::utils::bit_ops::test::BitsParam< 8>
//	 ,
//    ::bliss::utils::bit_ops::test::BitsParam<16> ,
//    ::bliss::utils::bit_ops::test::BitsParam<32>
> BitReverseLongArrayTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseLongArrayTest, BitReverseLongArrayTestTypes);





