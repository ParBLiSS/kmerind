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
 */



#include "utils/logging.h"

// include google test
#include <gtest/gtest.h>
#include "utils/bitgroup_ops.hpp"

#include <random>
#include <cstdint>

#include "common/kmer.hpp"
#include "common/test/kmer_reverse_helper.hpp"

#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "utils/timer.hpp"

// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseBenchmark : public ::testing::Test {
  protected:

    static T kmer;
    static T rev_gold;
    static T revcomp_gold;

    static bliss::common::test::KmerReverseHelper<T> helper;

    static constexpr size_t iterations = 100000000;

    static std::vector<long> chars;

  public:
    static void SetUpTestCase()
    {
      srand(23);
      for (size_t i = 0; i < T::nWords; ++i) {
        kmer.getData()[i] = static_cast<typename T::KmerWordType>(static_cast<long>(rand()) << 32) | static_cast<long>(rand());
      }

      chars.resize(1001);
      ::std::generate(chars.begin(), chars.end(), []() {return (static_cast<long>(rand()) << 32) | static_cast<long>(rand()); });

    }

};


template <typename T>
T KmerReverseBenchmark<T>::kmer;
template <typename T>
T KmerReverseBenchmark<T>::rev_gold;
template <typename T>
T KmerReverseBenchmark<T>::revcomp_gold;

template <typename T>
bliss::common::test::KmerReverseHelper<T> KmerReverseBenchmark<T>::helper;

template <typename T>
constexpr size_t KmerReverseBenchmark<T>::iterations;

template <typename T>
std::vector<long> KmerReverseBenchmark<T>::chars;


// to save some typing.  note that func is not a functor nor lambda function, so using macro is easier than as a templated function.
#define TEST_REV(name, func, kmertype, gold) do { \
    TypeParam km, rev, tmp; \
    km = KmerReverseBenchmark<kmertype>::kmer; \
    \
    TIMER_START(km); \
    \
    for (size_t i = 0; i < KmerReverseBenchmark<kmertype>::iterations; ++i) { \
      \
      tmp = func; \
      \
      rev ^= tmp; \
      \
      km.getData()[i % TypeParam::nWords] = static_cast<typename TypeParam::KmerWordType>(KmerReverseBenchmark<kmertype>::chars[i % 1001]) ; \
    } \
    TIMER_END(km, name, KmerReverseBenchmark<kmertype>::iterations); \
    \
    if (rev == km) \
     std::cout << "rev is same as km.  unlikely event." << std::endl; \
/*    if (rev != KmerReverseBenchmark<kmertype>::gold) \
      std::cout << "rev " << rev << std::endl << "gold " << KmerReverseBenchmark<kmertype>::gold << std::endl; \
      ASSERT_TRUE(rev == KmerReverseBenchmark<kmertype>::gold); \
*/    } while (0)


// to save some typing.  note that func is not a functor nor lambda function, so using macro is easier than as a templated function.
#define TEST_REV_BITOPS(name, simd, kmertype) do { \
    TypeParam km, rev, tmp; \
    km = KmerReverseBenchmark<kmertype>::kmer; \
    \
    TIMER_START(km); \
    \
    for (size_t i = 0; i < KmerReverseBenchmark<kmertype>::iterations; ++i) { \
      \
      memset(tmp.getData(), 0, TypeParam::nBytes); \
      \
      bliss::utils::bit_ops::reverse<TypeParam::bitsPerChar, simd>(tmp.getDataRef(), km.getDataRef()); \
      tmp.right_shift_bits(TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits);  \
      /* shift by remainder/padding. */ \
      \
      rev ^= tmp; \
      \
      km.getData()[i % TypeParam::nWords] = static_cast<typename TypeParam::KmerWordType>(KmerReverseBenchmark<kmertype>::chars[i % 1001]) ; \
    } \
    TIMER_END(km, name, KmerReverseBenchmark<kmertype>::iterations); \
    \
    if (rev == km) \
     std::cout << "rev is same as km.  unlikely event." << std::endl; \
    } while (0)


// to save some typing.  note that func is not a functor nor lambda function, so using macro is easier than as a templated function.
#define TEST_REVC_BITOPS(name, simd, kmertype) do { \
    TypeParam km, rev, tmp; \
    km = KmerReverseBenchmark<kmertype>::kmer; \
    \
    TIMER_START(km); \
    \
    for (size_t i = 0; i < KmerReverseBenchmark<kmertype>::iterations; ++i) { \
      \
      memset(tmp.getData(), 0, TypeParam::nBytes); \
      \
      switch (TypeParam::bitsPerChar) { \
        case 2: \
          bliss::utils::bit_ops::reverse<TypeParam::bitsPerChar, simd>(tmp.getDataRef(), km.getDataRef()); \
          bliss::utils::bit_ops::negate(tmp.getDataRef(), tmp.getDataRef()); \
          break; \
        case 3: \
        case 4: \
          bliss::utils::bit_ops::reverse<1, simd>(tmp.getDataRef(), km.getDataRef()); \
          break; \
        default: \
          break; \
      } \
      \
      tmp.right_shift_bits(TypeParam::nWords * sizeof(typename TypeParam::KmerWordType) * 8 - TypeParam::nBits);  \
      /* shift by remainder/padding. */ \
      \
      rev ^= tmp; \
      \
      km.getData()[i % TypeParam::nWords] = static_cast<typename TypeParam::KmerWordType>(KmerReverseBenchmark<kmertype>::chars[i % 1001]) ; \
    } \
    TIMER_END(km, name, KmerReverseBenchmark<kmertype>::iterations); \
    \
    if (rev == km) \
     std::cout << "rev is same as km.  unlikely event." << std::endl; \
    } while (0)



// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseBenchmark);

TYPED_TEST_P(KmerReverseBenchmark, reverse)
{
  TIMER_INIT(km);

  //TEST_REV("oldseq", KmerReverseBenchmark<TypeParam>::helper.reverse_serial(km), TypeParam, rev_gold);

  if (((TypeParam::bitsPerChar & (TypeParam::bitsPerChar - 1)) == 0) &&
      (::std::is_same<typename TypeParam::KmerAlphabet, bliss::common::DNA>::value ||
          ::std::is_same<typename TypeParam::KmerAlphabet, bliss::common::RNA>::value ||
           ::std::is_same<typename TypeParam::KmerAlphabet, bliss::common::DNA16>::value)) {
    TEST_REV("bswap", KmerReverseBenchmark<TypeParam>::helper.reverse_bswap(km), TypeParam, rev_gold);
    TEST_REV("swar", KmerReverseBenchmark<TypeParam>::helper.reverse_swar(km), TypeParam, rev_gold);

#ifdef __SSSE3__
    TEST_REV("ssse3", KmerReverseBenchmark<TypeParam>::helper.reverse_simd(km), TypeParam, rev_gold);
#endif

  }  // alphabet for DNA, RNA, and DNA16 are the only ones accelerated with simd type operations.

  TEST_REV_BITOPS("swar_new", ::bliss::utils::bit_ops::BIT_REV_SWAR, TypeParam);
#ifdef __SSSE3__
    TEST_REV_BITOPS("ssse3_new", ::bliss::utils::bit_ops::BIT_REV_SSSE3, TypeParam);
#endif
#ifdef __AVX2__
    TEST_REV_BITOPS("avx2_new", ::bliss::utils::bit_ops::BIT_REV_AVX2, TypeParam);
#endif
//  TEST_REV_BITOPS("seq_new", ::bliss::utils::bit_ops::BIT_REV_SEQ, TypeParam);
    TEST_REV("rev", km.reverse(), TypeParam, rev_gold);

  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}


TYPED_TEST_P(KmerReverseBenchmark, revcomp)
{
  TIMER_INIT(km);

  //TEST_REV("oldseqC", KmerReverseBenchmark<TypeParam>::helper.reverse_complement_serial(km), TypeParam, revcomp_gold);

  if (((TypeParam::bitsPerChar & (TypeParam::bitsPerChar - 1)) == 0) &&
      (::std::is_same<typename TypeParam::KmerAlphabet, bliss::common::DNA>::value ||
                ::std::is_same<typename TypeParam::KmerAlphabet, bliss::common::RNA>::value ||
                 ::std::is_same<typename TypeParam::KmerAlphabet, bliss::common::DNA16>::value) ) {
    TEST_REV("bswapC", KmerReverseBenchmark<TypeParam>::helper.reverse_complement_bswap(km), TypeParam, revcomp_gold);
    TEST_REV("swarC", KmerReverseBenchmark<TypeParam>::helper.reverse_complement_swar(km), TypeParam, revcomp_gold);

#ifdef __SSSE3__
    TEST_REV("ssse3C", KmerReverseBenchmark<TypeParam>::helper.reverse_complement_simd(km), TypeParam, revcomp_gold);
#endif

  }  // alphabet for DNA, RNA, and DNA16 are the only ones accelerated with simd type operations.


  TEST_REVC_BITOPS("swarC_new", ::bliss::utils::bit_ops::BIT_REV_SWAR, TypeParam);
#ifdef __SSSE3__
    TEST_REVC_BITOPS("ssse3_new", ::bliss::utils::bit_ops::BIT_REV_SSSE3, TypeParam);
#endif
#ifdef __AVX2__
    TEST_REVC_BITOPS("avx2_new", ::bliss::utils::bit_ops::BIT_REV_AVX2, TypeParam);
#endif
    // macro does not support bit_ops::BIT_REV_SEQ
    TEST_REV("revC", km.reverse_complement(), TypeParam, revcomp_gold);

  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}





//REGISTER_TYPED_TEST_CASE_P(KmerReverseBenchmark, rev_seq, rev_seq2, revcomp_seq, rev_bswap, revcomp_bswap, rev_swar, revcomp_swar, rev, revcomp, rev_ssse3, revcomp_ssse3);

REGISTER_TYPED_TEST_CASE_P(KmerReverseBenchmark, reverse, revcomp);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 15, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 47, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<192, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 15, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 47, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<192, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 47, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<192, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA16, uint64_t>
//    ::bliss::common::Kmer< 15, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer< 16, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer< 32, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer< 47, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer< 64, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer< 96, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer<128, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer<192, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer<256, bliss::common::DNA_IUPAC, uint64_t>,
//    ::bliss::common::Kmer< 7, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer< 15, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer< 16, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer< 32, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer< 47, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer< 64, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer< 96, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer<128, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer<192, bliss::common::ASCII, uint64_t>,
//    ::bliss::common::Kmer<256, bliss::common::ASCII, uint64_t>
> KmerReverseBenchmarkTypes;
//typedef ::testing::Types<
//    ::bliss::common::Kmer< 192, bliss::common::DNA,   uint64_t>,
//     ::bliss::common::Kmer< 96, bliss::common::DNA16,   uint64_t>
//> KmerReverseBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerReverseBenchmark, KmerReverseBenchmarkTypes);







// include files to test

//TESTS: Sequential, SWAR/BSWAP, SSSE3, AVX2 versions of kmer reverse.
//TESTS: for each, test kmer size, different word type, and Alphabet (bit group size)
//       different offsets, different bit group sizes, different word types, and different byte array lengths.

template <typename T>
class KmerReverseOpBenchmark : public ::testing::Test {
  protected:

    static T kmer;
    static T rev_gold;
    static T revcomp_gold;

    static bliss::common::test::KmerReverseHelper<T> helper;

    static constexpr size_t iterations = 100000000;

    static std::vector<long> chars;

  public:
    static void SetUpTestCase()
    {
      srand(23);
      for (size_t i = 0; i < T::nWords; ++i) {
        kmer.getData()[i] = static_cast<typename T::KmerWordType>(static_cast<long>(rand()) << 32) | static_cast<long>(rand());
      }

      chars.resize(1001);
      ::std::generate(chars.begin(), chars.end(), []() {return (static_cast<long>(rand()) << 32) | static_cast<long>(rand()); });

    }

    template <unsigned char SIMDType, typename TT = T>
    TT benchmark() {
      TT km, rev, tmp;
      km = kmer;


      for (size_t i = 0; i < iterations; ++i) {

        memset(tmp.getData(), 0, TT::nBytes);

        bliss::utils::bit_ops::reverse<TT::bitsPerChar, SIMDType>(tmp.getDataRef(), km.getDataRef());
        tmp.right_shift_bits(TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits);
        /* shift by remainder/padding. */

        rev ^= tmp;

        km.getData()[i % TT::nWords] = static_cast<typename TT::KmerWordType>(chars[i % 1001]) ;
      }

      return rev;
    }

    template <typename SIMDType, typename TT = T>
    TT benchmark() {
      TT km, rev, tmp;
      km = this->kmer;

      using MachWord = typename SIMDType::MachineWord;

      constexpr size_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;
      bliss::utils::bit_ops::bitgroup_ops<TT::bitsPerChar, SIMDType::SIMDVal> op;

      for (size_t i = 0; i < iterations; ++i) {

        bliss::utils::bit_ops::reverse<TT::bitsPerChar, SIMDType, shift>(tmp.getDataRef(), km.getDataRef(),
            [&op](MachWord const & src, MachWord & prev) {
          MachWord rev = op.reverse(src);
          MachWord res = bliss::utils::bit_ops::bit_or(prev, (bliss::utils::bit_ops::srli<(shift % 8)>(rev)));
          prev = rev; // bliss::utils::bit_ops::slli<(sizeof(MachWord) * 8 - shift)>(rev);
          return res;
        });

        rev ^= tmp;

        km.getData()[i % TT::nWords] = static_cast<typename TT::KmerWordType>(chars[i % 1001]) ;
      }

      return rev;
    }


    template <unsigned char SIMDType, typename TT = T>
    TT benchmark_c() {


      TT km, rev, tmp;
      km = KmerReverseOpBenchmark<TT>::kmer;

      for (size_t i = 0; i < iterations; ++i) {

        memset(tmp.getData(), 0, TT::nBytes);

        switch (TT::bitsPerChar) {
          case 2:
            bliss::utils::bit_ops::reverse<TT::bitsPerChar, SIMDType>(tmp.getDataRef(), km.getDataRef());
            bliss::utils::bit_ops::negate(tmp.getDataRef(), tmp.getDataRef());
            break;
          case 3:
          case 4:
            bliss::utils::bit_ops::reverse<1, SIMDType>(tmp.getDataRef(), km.getDataRef());
            break;
          default:
            break;
        }

        tmp.right_shift_bits(TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits);
        /* shift by remainder/padding. */

        rev ^= tmp;

        km.getData()[i % TT::nWords] = static_cast<typename TT::KmerWordType>(chars[i % 1001]) ;
      }

      return rev;
   }

    template <typename SIMDType, typename TT = T,
        typename ::std::enable_if<::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::RNA>::value, int>::type = 1>
    TT benchmark_c() {
      TT km, rev, tmp;
      km = KmerReverseOpBenchmark<TT>::kmer;
      using MachWord = typename SIMDType::MachineWord;

      constexpr size_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;
      bliss::utils::bit_ops::bitgroup_ops<TT::bitsPerChar, SIMDType::SIMDVal> op;

      for (size_t i = 0; i < KmerReverseOpBenchmark<TT>::iterations; ++i) {

        bliss::utils::bit_ops::reverse<TT::bitsPerChar, SIMDType, shift>(tmp.getDataRef(), km.getDataRef(),
            [&op](MachWord const & src, MachWord & prev) {
          MachWord rev = op.reverse(src);
          rev = bliss::utils::bit_ops::negate(rev);

          MachWord res = bliss::utils::bit_ops::bit_or(prev, (bliss::utils::bit_ops::srli<(shift % 8)>(rev)));

          prev = rev; //bliss::utils::bit_ops::slli<(sizeof(MachWord) * 8 - shift)>(rev);  // shift by full size returns unshifted data.

          return res;
        });

        rev ^= tmp;

        km.getData()[i % TT::nWords] = static_cast<typename TT::KmerWordType>(chars[i % 1001]) ;
      }

      return rev;
    }
    template <typename SIMDType, typename TT = T,
        typename ::std::enable_if<::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA16>::value ||
        ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::DNA6>::value ||
         ::std::is_same<typename TT::KmerAlphabet, ::bliss::common::RNA6>::value, int>::type = 1>
    TT benchmark_c() {
      TT km, rev, tmp;
      km = KmerReverseOpBenchmark<TT>::kmer;

      using MachineWord = typename SIMDType::MachineWord;

      constexpr size_t shift = TT::nWords * sizeof(typename TT::KmerWordType) * 8 - TT::nBits;
      bliss::utils::bit_ops::bitgroup_ops<1, SIMDType::SIMDVal> op1;

      for (size_t i = 0; i < KmerReverseOpBenchmark<TT>::iterations; ++i) {

        bliss::utils::bit_ops::reverse<1, SIMDType, shift>(tmp.getDataRef(), km.getDataRef(),
            [&op1](MachineWord const & src, MachineWord & prev) {
          MachineWord rev = op1.reverse(src);

          MachineWord res = bliss::utils::bit_ops::bit_or(prev, (bliss::utils::bit_ops::srli<(shift % 8)>(rev)));

          prev = rev; //bliss::utils::bit_ops::slli<(sizeof(MachineWord) * 8 - shift)>(rev);  // shift by full size returns unshifted data.

          return res;
        });

        rev ^= tmp;

        km.getData()[i % TT::nWords] = static_cast<typename TT::KmerWordType>(chars[i % 1001]) ;
      }

      return rev;
    }


};


template <typename T>
T KmerReverseOpBenchmark<T>::kmer;
template <typename T>
T KmerReverseOpBenchmark<T>::rev_gold;
template <typename T>
T KmerReverseOpBenchmark<T>::revcomp_gold;

template <typename T>
bliss::common::test::KmerReverseHelper<T> KmerReverseOpBenchmark<T>::helper;

template <typename T>
constexpr size_t KmerReverseOpBenchmark<T>::iterations;

template <typename T>
std::vector<long> KmerReverseOpBenchmark<T>::chars;
// indicate this is a typed test
TYPED_TEST_CASE_P(KmerReverseOpBenchmark);

TYPED_TEST_P(KmerReverseOpBenchmark, reverse_swar)
{
  TIMER_INIT(km);

  TIMER_START(km);
  TypeParam rev = this->template benchmark<bliss::utils::bit_ops::BIT_REV_SWAR>();
  TIMER_END(km, "rev swar", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "rev is same as kmer.  unlikely event." << std::endl;


  TIMER_START(km);
  rev = this->template benchmark<bliss::utils::bit_ops::BITREV_SWAR>();
  TIMER_END(km, "revop swar", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "rev is same as kmer.  unlikely event." << std::endl;


  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}

TYPED_TEST_P(KmerReverseOpBenchmark, reverse_ssse3)
{
  TIMER_INIT(km);

#ifdef __SSSE3__

  TIMER_START(km);
  TypeParam rev = this->template benchmark<bliss::utils::bit_ops::BIT_REV_SSSE3>();
  TIMER_END(km, "rev ssse3", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "rev is same as kmer.  unlikely event." << std::endl;


  TIMER_START(km);
  rev = this->template benchmark<bliss::utils::bit_ops::BITREV_SSSE3>();
  TIMER_END(km, "revop ssse3", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "rev is same as kmer.  unlikely event." << std::endl;

#endif

  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}
TYPED_TEST_P(KmerReverseOpBenchmark, reverse_avx2)
{
  TIMER_INIT(km);

#ifdef __AVX2__

  TIMER_START(km);
  TypeParam rev = this->template benchmark<bliss::utils::bit_ops::BIT_REV_AVX2>();
  TIMER_END(km, "rev avx", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "rev is same as kmer.  unlikely event." << std::endl;


  TIMER_START(km);
  rev = this->template benchmark<bliss::utils::bit_ops::BITREV_AVX2>();
  TIMER_END(km, "revop avx", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "rev is same as kmer.  unlikely event." << std::endl;

#endif


  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}


TYPED_TEST_P(KmerReverseOpBenchmark, revcomp_swar)
{
  TIMER_INIT(km);

  TIMER_START(km);
  TypeParam rev = this->template benchmark_c<bliss::utils::bit_ops::BIT_REV_SWAR>();
  TIMER_END(km, "revc swar", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "revc is same as kmer.  unlikely event." << std::endl;


  TIMER_START(km);
  rev = this->template benchmark_c<bliss::utils::bit_ops::BITREV_SWAR>();
  TIMER_END(km, "revopc swar", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "revc is same as kmer.  unlikely event." << std::endl;

  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}

TYPED_TEST_P(KmerReverseOpBenchmark, revcomp_ssse3)
{
  TIMER_INIT(km);


#ifdef __SSSE3__
  TIMER_START(km);
  TypeParam rev = this->template benchmark_c<bliss::utils::bit_ops::BIT_REV_SSSE3>();
  TIMER_END(km, "revc ssse3", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "revc is same as kmer.  unlikely event." << std::endl;


  TIMER_START(km);
  rev = this->template benchmark_c<bliss::utils::bit_ops::BITREV_SSSE3>();
  TIMER_END(km, "revopc ssse3", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "revc is same as kmer.  unlikely event." << std::endl;

#endif



  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}


TYPED_TEST_P(KmerReverseOpBenchmark, revcomp_avx2)
{
  TIMER_INIT(km);

#ifdef __AVX2__
  TIMER_START(km);
  TypeParam rev = this->template benchmark_c<bliss::utils::bit_ops::BIT_REV_AVX2>();
  TIMER_END(km, "revc avx", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "revc is same as kmer.  unlikely event." << std::endl;


  TIMER_START(km);
  rev = this->template benchmark_c<bliss::utils::bit_ops::BITREV_AVX2>();
  TIMER_END(km, "revopc avx", KmerReverseOpBenchmark<TypeParam>::iterations);

  if (rev == this->kmer)
   std::cout << "revc is same as kmer.  unlikely event." << std::endl;

#endif


  TIMER_REPORT(km, TypeParam::KmerAlphabet::SIZE);

}

REGISTER_TYPED_TEST_CASE_P(KmerReverseOpBenchmark, reverse_swar, reverse_ssse3, reverse_avx2, revcomp_swar, revcomp_ssse3, revcomp_avx2);

//////////////////// RUN the tests with different types.

// max of 50 cases
typedef ::testing::Types<
    ::bliss::common::Kmer< 15, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 47, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<192, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA,   uint64_t>,
    ::bliss::common::Kmer< 15, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 47, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<192, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA5,  uint64_t>,
    ::bliss::common::Kmer< 15, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 16, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 32, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 47, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 64, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer< 96, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<128, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<192, bliss::common::DNA16, uint64_t>,
    ::bliss::common::Kmer<256, bliss::common::DNA16, uint64_t>
> KmerReverseOpBenchmarkTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Bliss, KmerReverseOpBenchmark, KmerReverseOpBenchmarkTypes);

