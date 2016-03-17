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



template <typename T>
class BitReverseTest : public ::testing::Test {
  protected:

    BitReverseTestHelper<T::bitsPerGroup> helper;

    template <typename WORD_TYPE>
    bool is_reverse(WORD_TYPE const & orig, WORD_TYPE const & rev, uint8_t bit_offset = 0) {

      return helper.is_reverse(orig, rev, bit_offset);
    }

    bool is_reverse(uint8_t *out, uint8_t const * in, size_t len, uint8_t bit_offset = 0)  {

      return helper.is_reverse(out, in, len, bit_offset);
    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(BitReverseTest);

TYPED_TEST_P(BitReverseTest, reverse_uint8)
{
  using WORD_TYPE = uint8_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

      WORD_TYPE in;
  	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));


      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_uint16)
{
  using WORD_TYPE = uint16_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

        WORD_TYPE in;
    	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));

      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_uint32)
{
  using WORD_TYPE = uint32_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

        WORD_TYPE in;
    	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));

      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}


TYPED_TEST_P(BitReverseTest, reverse_uint64)
{
  using WORD_TYPE = uint64_t;

  TypeParam op;

  if (TypeParam::bitsPerGroup < (sizeof(WORD_TYPE) * 8) ) {
    for (size_t k = 0; k <= (32 - sizeof(WORD_TYPE)); ++k) {

        WORD_TYPE in;
    	  memcpy(&in, this->helper.input + k, sizeof(WORD_TYPE));

      WORD_TYPE out = op.reverse(in);

      ASSERT_TRUE(this->is_reverse(in, out));
    }
  }  // else too large, so don't do the test.

}

TYPED_TEST_P(BitReverseTest, reverse_short_array)
{

  TypeParam op;

  uint8_t BLISS_ALIGNED_ARRAY(out, 32, 32);

  unsigned int max = 8;

  for (unsigned int i = 1; i <= max; ++i ) {
    if (((i * 8) % TypeParam::bitsPerGroup) != 0) continue;  // i has to be a multiple of bytes for bitsPerGroup.


    if (TypeParam::bitsPerGroup == 3) {
      for (unsigned int k = 0; k <= (32 - i); ++k) {
        for (uint8_t j = 0; j < 8; ++j) {
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
      for (unsigned int k = 0; k <= (32 - i); ++k) {
        memset(out, 0, 32);

        op.reverse(out, this->helper.input + k, i);

        bool same = this->is_reverse(this->helper.input + k, out, i);

        if (!same) {
          std::cout << "in (MSB to LSB): ";
          for (int64_t j = i-1; j >= 0; --j) {
            std::cout << std::hex << static_cast<size_t>(this->helper.input[k + i - 1 - j]) << " ";

          }
          std::cout << std::endl;

          std::cout << "out (MSB to LSB): ";
          for (int64_t j = i-1; j >= 0; --j) {
            std::cout << std::hex << static_cast<size_t>(out[i - 1 - j]) << " ";

          }
          std::cout << std::endl;

          printf("array size = %d\n", i);
        }

        EXPECT_TRUE(same);
      }
    }
  }
}

// now register the test cases
REGISTER_TYPED_TEST_CASE_P(BitReverseTest,
                           reverse_uint8,
                           reverse_uint16,
                           reverse_uint32,
                           reverse_uint64,
                           reverse_short_array);


//////////////////// RUN the tests with different types.

typedef ::testing::Types<
    ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
     ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
      ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
       ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
        ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
         ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
          ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_SEQ>,
           ::bliss::utils::bit_ops::bitgroup_ops< 1, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
            ::bliss::utils::bit_ops::bitgroup_ops< 2, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
             ::bliss::utils::bit_ops::bitgroup_ops< 3, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
              ::bliss::utils::bit_ops::bitgroup_ops< 4, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
               ::bliss::utils::bit_ops::bitgroup_ops< 8, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
                ::bliss::utils::bit_ops::bitgroup_ops<16, ::bliss::utils::bit_ops::BIT_REV_SWAR> ,
                 ::bliss::utils::bit_ops::bitgroup_ops<32, ::bliss::utils::bit_ops::BIT_REV_SWAR>
> BitReverseTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, BitReverseTest, BitReverseTestTypes);


