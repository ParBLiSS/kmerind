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

/**
 * @file    test_kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 */



#include "utils/logging.h"

// include google test
#include <gtest/gtest.h>
#include "utils/bitgroup_ops.hpp"

#include <random>
#include <cstdint>

#include <atomic>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"

#include "common/test/kmer_reverse_helper.hpp"

// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseOpTest : public ::testing::Test {
  protected:

    T kmer;
    bliss::common::test::KmerReverseHelper<T> helper;

    static const size_t iterations = 100000;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {

        kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }
    }

    template <typename WORD_TYPE>
    static void print(WORD_TYPE const & word) {
      uint64_t const *ptr = reinterpret_cast<uint64_t const *>(&word);
      for (int i = (sizeof(WORD_TYPE) / sizeof(uint64_t)) - 1; i >= 0; --i) {
        std::cout << std::hex << std::setw(16) << std::setfill('0') << ptr[i] << " ";
      }
    }

    template <unsigned int BITS, unsigned char SIMD>
    struct reverse_op {
    	::bliss::utils::bit_ops::bitgroup_ops<BITS, SIMD> op;

    	template <typename WORD_TYPE>
    	inline WORD_TYPE operator()(WORD_TYPE const & src) const {
    		return op.reverse(src);
    	}
    };
    template <unsigned int BITS, unsigned char SIMD>
    struct reverse_negate_op {
    	::bliss::utils::bit_ops::bitgroup_ops<BITS, SIMD> op;

    	template <typename WORD_TYPE>
    	inline WORD_TYPE operator()(WORD_TYPE const & src) const {
    		return bliss::utils::bit_ops::bit_not(op.reverse(src));
    	}
    };
    template <unsigned int BITS, unsigned char SIMD>
    struct reverse_1bit_op {
    	::bliss::utils::bit_ops::bitgroup_ops<1, SIMD> op;

    	template <typename WORD_TYPE>
    	inline WORD_TYPE operator()(WORD_TYPE const & src) const {
    		return op.reverse(src);
    	}
    };

    template <typename SIMDType, typename TT = T,
        typename ::std::enable_if<::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA16>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::RNA>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA6>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::RNA6>::value, int>::type = 1>
    void test () {
      TT km, rev, rev_op;
      km = this->kmer;

      bool rev_same = true;
      bool local_rev_same = true;

      constexpr unsigned char SIMDVal = SIMDType::SIMDVal;

      constexpr uint16_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;
      constexpr uint16_t overlap = sizeof(typename SIMDType::MachineWord) % TT::bitsPerChar;

//      std::cout << "bit pow2 " << (uint64_t)shift << std::endl;
      reverse_op<TT::bitsPerChar, SIMDVal> op;


      for (size_t i = 0; i < this->iterations; ++i) {
        rev = km.reverse();

        bliss::utils::bit_ops::reverse_transform<SIMDType, shift, overlap,
          typename TT::KmerWordType, TT::nWords >(rev_op.getDataRef(), km.getDataRef(),
            op);

        local_rev_same = (rev == rev_op);

        if (!local_rev_same) {
          BL_WARNING("ERROR: rev diff at iter " << i << ":\n\tinput\t" << km << "\n\toutput\t" << rev << "\n\tgold\t" << rev_op);
        }
        rev_same &= local_rev_same;
        ASSERT_TRUE(local_rev_same);

        km.nextFromChar(rand() % TT::KmerAlphabet::SIZE);
      }
      EXPECT_TRUE(rev_same);

    }

    template <typename SIMDType, typename TT = T,
        typename ::std::enable_if<::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::RNA>::value, int>::type = 1>
    void testc () {
      TT km, rev, rev_op;
      km = this->kmer;

      bool rev_same = true;
      bool local_rev_same = true;

      constexpr unsigned char SIMDVal = SIMDType::SIMDVal;

      constexpr uint16_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;
      constexpr uint16_t overlap = sizeof(typename SIMDType::MachineWord) % TT::bitsPerChar;

//      std::cout << "bit pow2 " << (uint64_t)shift << std::endl;
      reverse_negate_op<TT::bitsPerChar, SIMDVal> op;


      for (size_t i = 0; i < this->iterations; ++i) {
        rev = km.reverse_complement();


        bliss::utils::bit_ops::reverse_transform<SIMDType, shift, overlap, typename TT::KmerWordType, TT::nWords>(rev_op.getDataRef(), km.getDataRef(),
            op);

        local_rev_same = (rev == rev_op);

        if (!local_rev_same) {
//          BL_WARNINGF("ERROR: rev diff at iter %lu:\n\tinput\t%s\n\toutput\t%s\n\tgold\t%s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_op.toAlphabetString().c_str());
            BL_WARNING("ERROR: rev diff at iter " << i << ":\n\tinput\t" << km << "\n\toutput\t" << rev << "\n\tgold\t" << rev_op);
        }
        rev_same &= local_rev_same;
        ASSERT_TRUE(local_rev_same);

        km.nextFromChar(rand() % TT::KmerAlphabet::SIZE);
      }
      EXPECT_TRUE(rev_same);
    }

    template <typename SIMDType, typename TT = T,
        typename ::std::enable_if<::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA16>::value ||
          ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA6>::value ||
          ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::RNA6>::value, int>::type = 0>
    void testc () {
      TT km, rev, rev_op;
      km = this->kmer;

      bool rev_same = true;
      bool local_rev_same = true;

      constexpr unsigned char SIMDVal = SIMDType::SIMDVal;

      constexpr uint16_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;

//      std::cout << "bit 3 " << (uint64_t)shift << std::endl;
      reverse_1bit_op<TT::bitsPerChar, SIMDVal> op1;

      for (size_t i = 0; i < this->iterations; ++i) {
        rev = km.reverse_complement();


        bliss::utils::bit_ops::reverse_transform<SIMDType, shift, 0, typename TT::KmerWordType, TT::nWords>(rev_op.getDataRef(), km.getDataRef(),
            op1);

        local_rev_same = (rev == rev_op);

        if (!local_rev_same) {
//          BL_WARNINGF("ERROR: rev diff at iter %lu:\n\tinput\t%s\n\toutput\t%s\n\tgold\t%s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_op.toAlphabetString().c_str());
            BL_WARNING("ERROR: rev diff at iter " << i << ":\n\tinput\t" << km << "\n\toutput\t" << rev << "\n\tgold\t" << rev_op);
        }
        rev_same &= local_rev_same;
        ASSERT_TRUE(local_rev_same);

        km.nextFromChar(rand() % TT::KmerAlphabet::SIZE);
      }
      EXPECT_TRUE(rev_same);

    }

};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseOpTest);

TYPED_TEST_P(KmerReverseOpTest, reverse_swar)
{
  this->template test<bliss::utils::bit_ops::BITREV_SWAR>();
}

TYPED_TEST_P(KmerReverseOpTest, reverse_ssse3)
{
#ifdef __SSSE3__
  this->template test<bliss::utils::bit_ops::BITREV_SSSE3>();
#else
  BL_WARNINGF("SSSE3 is not enabled or not available.");
#endif
}

TYPED_TEST_P(KmerReverseOpTest, reverse_avx2)
{
#ifdef __AVX2__
  this->template test<bliss::utils::bit_ops::BITREV_AVX2>();
#else
  BL_WARNINGF("AVX2 is not enabled or not available.");
#endif
}

TYPED_TEST_P(KmerReverseOpTest, reverse_auto)
{
	constexpr size_t bytes = sizeof(typename TypeParam::KmerWordType) * TypeParam::nWords;
  this->template test<bliss::utils::bit_ops::BITREV_AUTO_AGGRESSIVE<bytes> >();
}

TYPED_TEST_P(KmerReverseOpTest, revcomp_swar)
{
  this->template testc<bliss::utils::bit_ops::BITREV_SWAR>();
}

TYPED_TEST_P(KmerReverseOpTest, revcomp_ssse3)
{
#ifdef __SSSE3__
  this->template testc<bliss::utils::bit_ops::BITREV_SSSE3>();
#else
  BL_WARNINGF("SSSE3 is not enabled or not available.");
#endif
}

TYPED_TEST_P(KmerReverseOpTest, revcomp_avx2)
{
#ifdef __AVX2__
  this->template testc<bliss::utils::bit_ops::BITREV_AVX2>();
#else
  BL_WARNINGF("AVX2 is not enabled or not available.");
#endif
}

TYPED_TEST_P(KmerReverseOpTest, revcomp_auto)
{
	constexpr size_t bytes = sizeof(typename TypeParam::KmerWordType) * TypeParam::nWords;
  this->template testc<bliss::utils::bit_ops::BITREV_AUTO_AGGRESSIVE<bytes> >();
}



REGISTER_TYPED_TEST_CASE_P(KmerReverseOpTest,  reverse_swar, reverse_ssse3, reverse_avx2, reverse_auto, revcomp_swar, revcomp_ssse3, revcomp_avx2, revcomp_auto);

