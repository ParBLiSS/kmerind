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


// include files to test
#include "utils/test/bit_reverse_test_helper.hpp"


#ifdef __AVX2__
template <typename T>
class BitReverseAVX2Test : public ::testing::Test {
  protected:

    BitReverseTestHelper<T::bitsPerGroup> helper;

    bool is_reverse(__m256i const & orig, __m256i const & rev, uint8_t bit_offset = 0) {

      return helper.is_reverse(orig, rev, bit_offset);
    }

    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len, uint8_t bit_offset = 0)  {

      return helper.is_reverse(out, in, len, bit_offset);
    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseAVX2Test);

TYPED_TEST_P(BitReverseAVX2Test, reverse_m256i)
{
  TypeParam op;

  if (TypeParam::bitsPerGroup < 256 ) {
    __m256i in = _mm256_loadu_si256((__m256i*)(this->helper.input));

    __m256i out = op.reverse(in);

    ASSERT_TRUE(this->is_reverse(in, out));
  }  // else too large, so don't do the test.

}


TYPED_TEST_P(BitReverseAVX2Test, reverse_short_array)
{

  TypeParam op;

  uint8_t BLISS_ALIGNED_ARRAY(out, 32, 32);

  int max = 32;


  for (int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


    if (TypeParam::bitsPerGroup == 3) {
      for (int k = 0; k <= (32 - i); ++k) {
        for (int j = 0; j < 8; ++j) {
          memset(out, 0, 32);

          op.reverse(out, this->helper.input + k, i, j);

          bool same = this->is_reverse(this->helper.input + k, out, i, j);

          if (!same) {
            std::cout << "in: ";
            for (int l = 0; l < i; ++l) {
              std::cout << std::hex << static_cast<size_t>(this->helper.input[k + i - 1 - l]) << " ";

            }
            std::cout << std::endl;

            std::cout << "out: ";
            for (int l = 0; l < i; ++l) {
              std::cout << std::hex << static_cast<size_t>(out[i - 1 - l]) << " ";

            }
            std::cout << std::endl;


            printf("array size = %d, offset = %d, input byte offset = %d \n", i, j, k);
          }

          ASSERT_TRUE(same);
        }
      }
    } else {
      for (int k = 0; k <= (32 - i); ++k) {
        memset(out, 0, 32);

        op.reverse(out, this->helper.input + k, i, 0);

        bool same = this->is_reverse(this->helper.input + k, out, i, 0);

        if (!same) {
          printf("array size = %d\n", i);
        }

        ASSERT_TRUE(same);
      }
    }
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseAVX2Test, reverse_m256i, reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
     ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
      ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
       ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
        ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
         ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
          ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
           ::bliss::utils::bit_ops::bitgroup_ops<64, ::bliss::utils::bit_ops::BIT_REV_AVX2> ,
            ::bliss::utils::bit_ops::bitgroup_ops<128, ::bliss::utils::bit_ops::BIT_REV_AVX2>
> BitReverseAVX2TestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseAVX2Test, BitReverseAVX2TestTypes);

#endif


