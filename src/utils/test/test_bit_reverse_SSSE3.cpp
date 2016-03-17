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

  uint8_t BLISS_ALIGNED_ARRAY(out, 32, 32);

  int max = 16;

  for (int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


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


