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
class KmerReverseTest : public ::testing::Test {
  protected:

    T kmer;
    bliss::common::test::KmerReverseHelper<T> helper;

    static const size_t iterations = 100000;

    virtual void SetUp()
    {
      srand(0);
      for (unsigned int i = 0; i < T::size; ++i) {
//    	  printf(" setup next from char\n");
        kmer.nextFromChar(rand() % T::KmerAlphabet::SIZE);
      }

    }

	template <typename data_type, size_t data_size>
    void print(data_type (&w)[data_size]) {
      std::cout << "data type size " << std::dec << sizeof(data_type) << " len " << data_size << ": ";
      for (int k = data_size -1 ; k >= 0; --k) {
        std::cout << std::hex << std::setfill('0') << std::setw(sizeof(data_type) * 2) << static_cast<size_t>(w[k]) << " ";
      }
      std::cout << std::endl;
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
      std::cout << "output: ";  this->print(km.getDataRef());
      std::cout << "gold: ";  this->print(gold.getDataRef());
    }
    ASSERT_TRUE(local_same);

    same &= local_same;

    gold.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
    km = gold;
  }
  EXPECT_TRUE(same);
}



TYPED_TEST_P(KmerReverseTest, reverse_seq)
{
  TypeParam km, rev, rev_seq, rev2;
  km = this->kmer;

  bool rev_same = true;
  bool local_rev_same = true;

  uint8_t* out = reinterpret_cast<uint8_t*>(rev.getDataRef());
  const uint8_t* in = reinterpret_cast<uint8_t const *>(km.getData());

  for (size_t i = 0; i < this->iterations; ++i) {
//	printf("rev_seq serial %lu\n", i);
    rev_seq = this->helper.reverse_serial(km);

//	printf("rev_seq BITSEQ reverse %lu\n", i);
    bliss::utils::bit_ops::reverse<TypeParam::bitsPerChar, bliss::utils::bit_ops::BIT_REV_SEQ>(out, in, TypeParam::nBytes);
//    bliss::utils::bit_ops::right_shift<bliss::utils::bit_ops::BITREV_AUTO_AGGRESSIVE<TypeParam::nWords * sizeof(typename TypeParam::KmerWordType)>,
//    		>(rev2.getDataRef(), rev.getDataRef());
//	printf("rev_seq shift %lu\n", i);

    ::std::atomic_thread_fence(::std::memory_order_seq_cst);
    rev.template right_shift_bits<(TypeParam::nBytes * 8 - TypeParam::nBits)>();  // shift by remainder bits..

    local_rev_same = (rev == rev_seq);

    if (!local_rev_same) {
    	printf("nwords %u word size %lu, shift by %lu\n",
    			TypeParam::nWords,
    			sizeof(typename TypeParam::KmerWordType),
    			(TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits));
      BL_ERROR("ERROR: rev diff at iter " << i << ":");
//      printf("printing km\n");
      BL_ERROR("\tinput " << km.toAlphabetString().c_str());
//      printf("printing rev2\n");
      BL_ERROR("\toutput " << rev.toAlphabetString().c_str());
//      printf("printing rev\n");
      BL_ERROR("\tgold " << rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    ASSERT_TRUE(local_rev_same);

    rev_same &= local_rev_same;

//	printf("rev_seq BITSEQ next from char %lu\n", i);
    km.nextFromChar(rand() % TypeParam::KmerAlphabet::SIZE);
  }
  ASSERT_TRUE(rev_same);

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
//		printf("rev_bswap rev-seq %lu\n", i);
    rev_seq = this->helper.reverse_serial(km);
//	printf("rev_bswap revcomp-seq %lu\n", i);
    revcomp_seq = this->helper.reverse_complement_serial(km);

//	printf("rev_bswap rev-bswap %lu\n", i);
    rev = this->helper.reverse_bswap(km);
//	printf("rev_bswap revcomp-bswap %lu\n", i);
    revcomp = this->helper.reverse_complement_bswap(km);

    local_rev_same = (rev == rev_seq);
    local_revcomp_same = (revcomp == revcomp_seq);

    if (!local_rev_same) {
//    	printf("rev_bswap printing %lu\n", i);

      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i,
    		  km.toAlphabetString().c_str(),
    		  rev.toAlphabetString().c_str(),
    		  rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
//      printf("rev_bswap printing done %lu\n", i);
    }
    if (!local_revcomp_same) {
//    	printf("rev_bswap printing %lu\n", i);
      BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i,
    		  km.toAlphabetString().c_str(),
    		  revcomp.toAlphabetString().c_str(),
    		  revcomp_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(revcomp.getDataRef());
      std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
//  	printf("rev_bswap printing done %lu\n", i);
    }
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

    rev_same &= local_rev_same;
    revcomp_same &= local_revcomp_same;

//	printf("rev_bswap next char %lu\n", i);
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
      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    if (!local_revcomp_same) {
        BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
        std::cout << "output: ";  this->print(revcomp.getDataRef());
        std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
    }
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

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
      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    if (!local_revcomp_same) {
        BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
        std::cout << "output: ";  this->print(revcomp.getDataRef());
        std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
    }
    ASSERT_TRUE(local_rev_same);
    ASSERT_TRUE(local_revcomp_same);

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
      BL_ERRORF("ERROR: rev diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), rev.toAlphabetString().c_str(), rev_seq.toAlphabetString().c_str());
      std::cout << "output: ";  this->print(rev.getDataRef());
      std::cout << "gold: ";  this->print(rev_seq.getDataRef());
    }
    if (!local_revcomp_same) {
        BL_ERRORF("ERROR: revcomp diff at iter %lu:\n\tinput %s\n\toutput %s\n\tgold %s", i, km.toAlphabetString().c_str(), revcomp.toAlphabetString().c_str(), revcomp_seq.toAlphabetString().c_str());
        std::cout << "output: ";  this->print(revcomp.getDataRef());
        std::cout << "gold: ";  this->print(revcomp_seq.getDataRef());
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


REGISTER_TYPED_TEST_CASE_P(KmerReverseOpTest, reverse_swar, reverse_ssse3, reverse_avx2, reverse_auto, revcomp_swar, revcomp_ssse3, revcomp_avx2, revcomp_auto);

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
