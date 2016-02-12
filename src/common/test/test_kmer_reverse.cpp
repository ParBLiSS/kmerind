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

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"

#include "common/test/kmer_reverse_helper.hpp"



// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseTest : public ::testing::Test {
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

};


// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseTest);

TYPED_TEST_P(KmerReverseTest, reverse_seq_self)
{
  TypeParam km, gold;
  km = this->kmer;
  gold = this->kmer;

  bool same = true;
  bool local_same;

  for (size_t i = 0; i < this->iterations; ++i) {
    km = this->helper.reverse_serial(km);
    km = this->helper.reverse_complement_serial(km);

    km = this->helper.reverse_serial(km);
    km = this->helper.reverse_complement_serial(km);

    local_same = (km == gold);

    if (!local_same) {
      BL_DEBUGF("ERROR: seq rev-revcomp-rev-revcomp diff at iter %lu:\n\tinput %s\n\toutput %s", i, gold.toAlphabetString().c_str(), km.toAlphabetString().c_str());
    }
    same &= local_same;

    gold.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
    km = gold;
  }
  EXPECT_TRUE(same);
}



TYPED_TEST_P(KmerReverseTest, reverse_seq)
{
  TypeParam km, rev, rev_seq;
  km = this->kmer;

  bool rev_same = true;
  bool local_rev_same = true;

  uint8_t* out = reinterpret_cast<uint8_t*>(rev.getData());
  const uint8_t* in = reinterpret_cast<uint8_t const *>(km.getData());

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);

    bliss::utils::bit_ops::reverse<TypeParam::bitsPerChar, bliss::utils::bit_ops::BIT_REV_SEQ>(out, in, TypeParam::nBytes);
    rev.right_shift_bits(TypeParam::nBytes * 8 - TypeParam::nBits);  // shift by remainder bits..

    local_rev_same = (rev == rev_seq);

    if (!local_rev_same) {
      BL_ERROR("ERROR: rev diff at iter " << i << ":\n\tinput " << km.toAlphabetString().c_str() << "\n\toutput " << rev.toAlphabetString().c_str() << "\n\tgold " << rev_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);

}

TYPED_TEST_P(KmerReverseTest, reverse_bswap)
{
  if (TypeParam::bitsPerChar == 3) return;

  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = this->helper.reverse_bswap(km);
    revcomp = this->helper.reverse_complement_bswap(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      BL_DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}

TYPED_TEST_P(KmerReverseTest, reverse_swar)
{
  if (TypeParam::bitsPerChar == 3) return;

  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = this->helper.reverse_swar(km);
    revcomp = this->helper.reverse_complement_swar(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      BL_DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}


TYPED_TEST_P(KmerReverseTest, reverse_ssse3)
{
#ifdef __SSSE3__
  if (TypeParam::bitsPerChar == 3) return;


  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = this->helper.reverse_simd(km);
    revcomp = this->helper.reverse_complement_simd(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      BL_DEBUGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;


    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);
#else
  BL_WARNINGF("SSSE3 is not enabled or not available.");
#endif

}

TYPED_TEST_P(KmerReverseTest, reverse)
{
  TypeParam km, rev, revcomp, rev_seq, revcomp_seq;
  km = this->kmer;

  bool rev_same = true;
  bool revcomp_same = true;
  bool local_rev_same = true;
  bool local_revcomp_same = true;

  for (size_t i = 0; i < this->iterations; ++i) {
    rev_seq = this->helper.reverse_serial(km);
    revcomp_seq = this->helper.reverse_complement_serial(km);

    rev = km.reverse();
    revcomp = km.reverse_complement();

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
      BL_DEBUGF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
    }
    if (!local_revcomp_same) {
      BL_WARNINGF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
    }
    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  EXPECT_TRUE(rev_same);
  EXPECT_TRUE(revcomp_same);

}




REGISTER_TYPED_TEST_CASE_P(KmerReverseTest, reverse_seq_self, reverse_seq, reverse_bswap, reverse_swar, reverse_ssse3, reverse);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA,   uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint32_t>,  // 1 word, full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint32_t>,  // 2 words, full
    ::bliss::common::Kmer< 40, bliss::common::DNA,   uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA,   uint32_t>,  // 3 words, full
    ::bliss::common::Kmer<  7, bliss::common::DNA,   uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA,   uint16_t>,  // 1 word, full
    ::bliss::common::Kmer<  9, bliss::common::DNA,   uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint16_t>,  // 2 words, full
    ::bliss::common::Kmer<  3, bliss::common::DNA,    uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA,    uint8_t>,  // 1 word, full
    ::bliss::common::Kmer<  5, bliss::common::DNA,    uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA,    uint8_t>,  // 2 words, full
    ::bliss::common::Kmer< 21, bliss::common::DNA5,  uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 22, bliss::common::DNA5,  uint64_t>,  // 2 word, not full
    ::bliss::common::Kmer< 42, bliss::common::DNA5,  uint64_t>,  // 2 words, not full
    ::bliss::common::Kmer< 43, bliss::common::DNA5,  uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,  // 3 words, full
    ::bliss::common::Kmer<  2, bliss::common::DNA5,   uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  3, bliss::common::DNA5,   uint8_t>,  // 2 word, not full
    ::bliss::common::Kmer<  5, bliss::common::DNA5,   uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  6, bliss::common::DNA5,   uint8_t>,  // 3 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA5,   uint8_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 40, bliss::common::DNA16, uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA16, uint64_t>,  // 3 words, full
    ::bliss::common::Kmer<  7, bliss::common::DNA16, uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA16, uint32_t>,  // 1 word, full
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint32_t>,  // 2 words, full
    ::bliss::common::Kmer< 20, bliss::common::DNA16, uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 24, bliss::common::DNA16, uint32_t>,  // 3 words, full
    ::bliss::common::Kmer<  3, bliss::common::DNA16, uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA16, uint16_t>,  // 1 word, full
    ::bliss::common::Kmer<  5, bliss::common::DNA16, uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA16, uint16_t>,  // 2 words, full
    ::bliss::common::Kmer<  1, bliss::common::DNA16,  uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  2, bliss::common::DNA16,  uint8_t>,  // 1 word, full
    ::bliss::common::Kmer<  3, bliss::common::DNA16,  uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA16,  uint8_t>,   // 2 words, full
//    ::bliss::common::Kmer< 15, bliss::common::DNA_IUPAC,   uint64_t>,  // 1 word, not full.  not tested since it is not a bijection.
//    ::bliss::common::Kmer< 40, bliss::common::DNA_IUPAC,   uint64_t>,  // 3 word, not full
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 21, bliss::common::ASCII,  uint8_t>  // 3 words, not full
> KmerReverseTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerReverseTest, KmerReverseTestTypes);



















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

      using MachWord = typename SIMDType::MachineWord;
      constexpr unsigned char SIMDVal = SIMDType::SIMDVal;

      constexpr unsigned char shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;

//      std::cout << "bit pow2 " << (uint64_t)shift << std::endl;
      bliss::utils::bit_ops::bitgroup_ops<TT::bitsPerChar, SIMDVal> op;


      for (size_t i = 0; i < this->iterations; ++i) {
        rev = km.reverse();

        bliss::utils::bit_ops::reverse<TT::bitsPerChar, SIMDType, shift>(rev_op.getDataRef(), km.getDataRef(),
            [&op](MachWord const & src) {
          return op.reverse(src);
        });

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

      using MachWord = typename SIMDType::MachineWord;
      constexpr unsigned char SIMDVal = SIMDType::SIMDVal;

      constexpr unsigned char shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;

//      std::cout << "bit pow2 " << (uint64_t)shift << std::endl;
      bliss::utils::bit_ops::bitgroup_ops<TT::bitsPerChar, SIMDVal> op;


      for (size_t i = 0; i < this->iterations; ++i) {
        rev = km.reverse_complement();


        bliss::utils::bit_ops::reverse<TT::bitsPerChar, SIMDType, shift>(rev_op.getDataRef(), km.getDataRef(),
            [&op](MachWord const & src) {
          return bliss::utils::bit_ops::negate(op.reverse(src));
        });

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

      using MachWord = typename SIMDType::MachineWord;
      constexpr unsigned char SIMDVal = SIMDType::SIMDVal;

      constexpr uint16_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;

//      std::cout << "bit 3 " << (uint64_t)shift << std::endl;
      bliss::utils::bit_ops::bitgroup_ops<1, SIMDVal> op1;

      for (size_t i = 0; i < this->iterations; ++i) {
        rev = km.reverse_complement();


        bliss::utils::bit_ops::reverse<1, SIMDType, shift>(rev_op.getDataRef(), km.getDataRef(),
            [&op1](MachWord const & src) {
          return op1.reverse(src);

        });

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


REGISTER_TYPED_TEST_CASE_P(KmerReverseOpTest, reverse_swar, reverse_ssse3, reverse_avx2, revcomp_swar, revcomp_ssse3, revcomp_avx2);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 31, bliss::common::DNA,   uint64_t>,  // 1 word, not full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,  // 1 word, full
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,  // 2 words, full
    ::bliss::common::Kmer< 80, bliss::common::DNA,   uint64_t>,  // 3 words, not full
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,  // 3 words, full
    ::bliss::common::Kmer< 15, bliss::common::DNA,   uint32_t>,  // 1 word, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint32_t>,  // 1 word, full
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint32_t>,  // 2 words, full
    ::bliss::common::Kmer< 40, bliss::common::DNA,   uint32_t>,  // 3 words, not full
    ::bliss::common::Kmer< 48, bliss::common::DNA,   uint32_t>,  // 3 words, full
    ::bliss::common::Kmer<  7, bliss::common::DNA,   uint16_t>,  // 1 word, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA,   uint16_t>,  // 1 word, full
    ::bliss::common::Kmer<  9, bliss::common::DNA,   uint16_t>,  // 2 words, not full
    ::bliss::common::Kmer< 16, bliss::common::DNA,   uint16_t>,  // 2 words, full
    ::bliss::common::Kmer<  3, bliss::common::DNA,    uint8_t>,  // 1 word, not full
    ::bliss::common::Kmer<  4, bliss::common::DNA,    uint8_t>,  // 1 word, full
    ::bliss::common::Kmer<  5, bliss::common::DNA,    uint8_t>,  // 2 words, not full
    ::bliss::common::Kmer<  8, bliss::common::DNA,    uint8_t>,  // 2 words, full
     ::bliss::common::Kmer< 21, bliss::common::DNA5,  uint64_t>,  // 1 word, not full
     ::bliss::common::Kmer< 22, bliss::common::DNA5,  uint64_t>,  // 2 word, not full
     ::bliss::common::Kmer< 42, bliss::common::DNA5,  uint64_t>,  // 2 words, not full
     ::bliss::common::Kmer< 43, bliss::common::DNA5,  uint64_t>,  // 3 words, not full
     ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,  // 3 words, full
     ::bliss::common::Kmer<  2, bliss::common::DNA5,   uint8_t>,  // 1 word, not full
     ::bliss::common::Kmer<  3, bliss::common::DNA5,   uint8_t>,  // 2 word, not full
     ::bliss::common::Kmer<  5, bliss::common::DNA5,   uint8_t>,  // 2 words, not full
     ::bliss::common::Kmer<  6, bliss::common::DNA5,   uint8_t>,  // 3 words, not full
     ::bliss::common::Kmer<  8, bliss::common::DNA5,   uint8_t>,  // 3 words, full
     ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>,  // 1 word, not full
     ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>,  // 1 word, full
     ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,  // 2 words, full
     ::bliss::common::Kmer< 40, bliss::common::DNA16, uint64_t>,  // 3 words, not full
     ::bliss::common::Kmer< 48, bliss::common::DNA16, uint64_t>,  // 3 words, full
     ::bliss::common::Kmer<  7, bliss::common::DNA16, uint32_t>,  // 1 word, not full
     ::bliss::common::Kmer<  8, bliss::common::DNA16, uint32_t>,  // 1 word, full
     ::bliss::common::Kmer< 16, bliss::common::DNA16, uint32_t>,  // 2 words, full
     ::bliss::common::Kmer< 20, bliss::common::DNA16, uint32_t>,  // 3 words, not full
     ::bliss::common::Kmer< 24, bliss::common::DNA16, uint32_t>,  // 3 words, full
     ::bliss::common::Kmer<  3, bliss::common::DNA16, uint16_t>,  // 1 word, not full
     ::bliss::common::Kmer<  4, bliss::common::DNA16, uint16_t>,  // 1 word, full
     ::bliss::common::Kmer<  5, bliss::common::DNA16, uint16_t>,  // 2 words, not full
     ::bliss::common::Kmer<  8, bliss::common::DNA16, uint16_t>,  // 2 words, full
     ::bliss::common::Kmer<  1, bliss::common::DNA16,  uint8_t>,  // 1 word, not full
     ::bliss::common::Kmer<  2, bliss::common::DNA16,  uint8_t>,  // 1 word, full
     ::bliss::common::Kmer<  3, bliss::common::DNA16,  uint8_t>,  // 2 words, not full
     ::bliss::common::Kmer<  4, bliss::common::DNA16,  uint8_t>   // 2 words, full

> KmerReverseOpTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerReverseOpTest, KmerReverseOpTestTypes);
