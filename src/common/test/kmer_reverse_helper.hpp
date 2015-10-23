/**
 * @file    kmer_reverse_helper.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the Kmer data type.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_TEST_KMER_REV_HELPER_H
#define BLISS_COMMON_TEST_KMER_REV_HELPER_H

// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <iterator>
#include <type_traits>
#include <typeinfo>   // typeid
#include <string>
#include <iostream>   //std::cout, for debugging.
#include <sstream>
#include <algorithm>
#include <utility>  // std::pair
#include <stack>
#include <cstring>  // memset, memcpy
#include <iomanip>  // setfill, setw, for std::cout


// own includes
#include "common/base_types.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "common/bit_ops.hpp"
#include "common/padding.hpp"
#include "utils/kmer_utils.hpp"
#include "utils/bitgroup_ops.hpp"


// #include <xmmintrin.h>  // sse2
#if defined(__SSSE3__)
#include <tmmintrin.h>  // ssse3  includes sse2
#endif

#include <bits/byteswap.h>


#define INLINE inline

namespace bliss
{

  namespace common
  {

    namespace test
    {

  /**
   * @brief   Implements a general templated k-mer class.
   *
   * The k-mer size (number of characters) is a template parameter, thus
   * the k-mer size for this k-mer needs to be fixed at compile time. The given
   * k-mer size and the `BITS_PER_CHAR` determine the size of the underlying
   * data array at compile time.
   *
   * A kmer's Most Significant Bit (MSB) corresponds to the prefix of the input string.
   *
   * @details  In memory organization of the kmer is described below.
   *      Kmer is packed bitwise based on the alphabet.  i.e. alphabet defines required
   *      number of bits per character.  characters in a kmer are packed together,
   *      with padding occurring at the most significant bits.  In case of multi-word
   *      kmer, array element 0 has LSB for the kmer, and the highest element has the MSB,
   *      as well as the padding.  The first character encounter in the sequence is more significant,
   *      so as to allow sorting by prefix.   For example:
   *        string:  ATCGGACTTA
   *        2 bits per DNA character, then we need 3 bytes, with 4 bits padding.
   *
   *        byte 2 in array:  00AT
   *        byte 1 in array:  CGGA
   *        byte 0 in array:  CTTA
   *
   *      shifting the kmer window to the next character in sequence results in overall
   *      left shift of bits in the array
   *
   *        byte 2 in array:  00TC
   *        byte 1 in array:  GGAC
   *        byte 0 in array:  TTAX
   *
   *      lexicographic comparison then starts with highest order
   *      prefix at MSB side.
   *
   *      Can generate the reverse kmer from unpacked characters. (complement is done outside of the kmer.  rationale is
   *      that it is more semantically clear)  Also, no reverse from packed characters, since the complement
   *      operation would be complicated on packed string and probably will not be used often anyways.
   *
   * @todo: implement and refer to the dynamic k-mer (which will be more
   *        inefficient)
   *
   * @tparam  KMER_SIZE   The size (number of characters) of the k-mer.
   * @tparam  ALPHABET    The Alphabet from which the characters in the k-mer are drawn.
   *                          E.g. DNA
   * @tparam  WORD_TYPE       The unsigned integer type to be used as the
   *                          storage type. This can be of any size from
   *                          uint8_t up to uint64_t. (default is uint64_t)
   */
  template <typename Kmer>
  class KmerReverseHelper
  {

   protected:
    using WORD_TYPE = typename Kmer::KmerWordType;
    unsigned int size = Kmer::size;
    using ALPHABET = typename Kmer::KmerAlphabet;


    using MACH_WORD_TYPE = size_t;

    static constexpr int stride = sizeof(MACH_WORD_TYPE) / sizeof(WORD_TYPE);
    static constexpr int iters = Kmer::nWords / stride;
    static constexpr int rem = Kmer::nWords % stride;

    static constexpr int simd_stride = 16 / sizeof(WORD_TYPE);
    static constexpr int simd_iters = Kmer::nWords / simd_stride;
    static constexpr int simd_rem = Kmer::nWords % simd_stride;


    static constexpr uint8_t simd_mask_b[16] alignas(16) = {0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F};

    // mask for reversing the order of items
    static constexpr uint8_t simd_rev_mask_b[16] alignas(16) = {0x0F,0x0E,0x0D,0x0C,0x0B,0x0A,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,0x00};

    // lookup table for reversing bits in a byte in groups of 1
    static constexpr uint8_t simd_lut1_lo_b[16] alignas(16) = {0x00,0x08,0x04,0x0c,0x02,0x0a,0x06,0x0e,0x01,0x09,0x05,0x0d,0x03,0x0b,0x07,0x0f};

    static constexpr uint8_t simd_lut1_hi_b[16] alignas(16) = {0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0};

    // lookup table for reversing bits in a byte in groups of 2
    static constexpr uint8_t simd_lut2_lo_b[16] alignas(16) = {0x00,0x04,0x08,0x0c,0x01,0x05,0x09,0x0d,0x02,0x06,0x0a,0x0e,0x03,0x07,0x0b,0x0f};

    static constexpr uint8_t simd_lut2_hi_b[16] alignas(16) = {0x00,0x40,0x80,0xc0,0x10,0x50,0x90,0xd0,0x20,0x60,0xa0,0xe0,0x30,0x70,0xb0,0xf0};

    // lookup table for reversing bits in a byte in groups of 4 - no need, just use the simd_mask.

  public:
  
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    INLINE Kmer reverse_serial(Kmer const & src) const
    {
      // TODO implement logarithmic version (logarithmic in number of bits in a word, linear in number of words).  non power of 2 bitsPerChar is tricky because of byte boundaries.
  
      /* Linear (inefficient) reverse: */
  
      // get temporary copy of this
      Kmer tmp_copy = src;
      Kmer result;
  
      // get lower most bits from the temp copy and push them into the lower bits
      // of this
      for (unsigned int i = 0; i < size; ++i)
      {
        result <<= 1;   // shifting the whole thing, inefficient but correct,
                                            // especially for char that cross word boundaries.
        // copy `bitsperChar` least significant bits
        copyBitsFixed<WORD_TYPE, Kmer::bitsPerChar>(result.getData()[0], tmp_copy.getConstData()[0]);

        tmp_copy >>= 1;
      }
  
      // result already was 0 to begin with, so no need to sanitize
//      // set ununsed bits to 0
//      this->do_sanitize();

      return result;
    }
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    INLINE Kmer reverse_complement_serial(Kmer const & src) const
    {
      // TODO implement logarithmic version (logarithmic in number of bits in a word, linear in number of words).  non power of 2 bitsPerChar is tricky because of byte boundaries.

      /* Linear (inefficient) reverse: */

      // get temporary copy of this
      Kmer tmp_copy = src;
      Kmer result;
      // get lower most bits from the temp copy and push them into the lower bits
      // of this
      WORD_TYPE tmp;

      for (unsigned int i = 0; i < size; ++i)
      {
        result <<= 1;   // shifting the whole thing, inefficient but correct,
                                            // especially for char that cross word boundaries.

        tmp = ALPHABET::TO_COMPLEMENT[tmp_copy.getConstData()[0] & getLeastSignificantBitsMask<WORD_TYPE>(Kmer::bitsPerChar)];

        // copy `bitsperChar` least significant bits
        copyBitsFixed<WORD_TYPE, Kmer::bitsPerChar>(result.getData()[0], tmp);
        tmp_copy >>= 1;
      }

      // result already was 0 to begin with, so no need to sanitize
//      // set ununsed bits to 0
//      this->do_sanitize();

      return result;

    }

#if defined(__SSSE3__)
    static void print_simd_register(__m128i const & v, std::string const & label) {
      WORD_TYPE tmp[simd_stride];
      _mm_storeu_si128((__m128i*)tmp, v);
      std::cout << label << ": ";
      size_t width = sizeof(WORD_TYPE) * 2;
      for (int i = simd_stride - 1; i >= 0; --i) {
        std::cout << std::hex << std::setfill('0') << std::setw(width) << tmp[i] << " ";
      }
      std::cout << std::dec << std::endl;
    }


    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.  8 to 100 times faster.
    // this code is inspired by http://stackoverflow.com/questions/746171/best-algorithm-for-bit-reversal-from-msb-lsb-to-lsb-msb-in-c
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
    ::std::is_same<A, bliss::common::RNA>::value, __m128i>::type word_reverse_simd(__m128i const & b) {

//        print(b, "b");

      // load constants:
      __m128i mask_lo = _mm_load_si128((__m128i*)simd_mask_b);

//        print(mask_lo, "mask_lo");

      //== first shuffle (reverse) the bytes
      __m128i r = _mm_shuffle_epi8(b, _mm_load_si128((__m128i*)simd_rev_mask_b));// *(reinterpret_cast<const __m128i*>(simd_rev_mask_b)));                                    // SSSE3

//        print(r, "r");

      //== get the lut indices.
      // lower 4 bits
      __m128i lo = _mm_and_si128(mask_lo, r); //*(reinterpret_cast<const __m128i*>(simd_mask_b)), r);                                          // SSE2
      // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
      __m128i hi = _mm_srli_epi16(_mm_andnot_si128(mask_lo, r), 4);   //*(reinterpret_cast<const __m128i*>(simd_mask_b)), r), 4);                    // SSE2

//        print(lo, "lo");
//        print(hi, "hi");

        __m128i slo = _mm_shuffle_epi8(_mm_load_si128((__m128i*)simd_lut2_hi_b), lo);
        __m128i shi = _mm_shuffle_epi8(_mm_load_si128((__m128i*)simd_lut2_lo_b), hi);

//        print(slo, "slo");
//        print(shi, "shi");

        //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
      __m128i res =  _mm_or_si128(slo, shi);                           // SSSE3, SSE2

//      print(res, "result");

      return res;

    }
    /// do reverse complement.  8 to 100 times faster than serially reversing the bits.
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
        ::std::is_same<A, bliss::common::RNA>::value, __m128i>::type word_reverse_complement_simd(__m128i const & b) {

      // DNA type, requires that alphabet be setup so that negation produces the complement.  cmpeq to return a value with all bits set to 1
      return _mm_xor_si128(word_reverse_simd<A>(b), _mm_cmpeq_epi8(b, b));  // sse does not have negation. use xor with FFFF   // SSE2
    }

    // reverse via bit swapping in 4 bit increment.  this complement table is then used for lookup.
    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.
    /// 8 to 50 times faster than sequentially reverse the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value, __m128i>::type word_reverse_simd(__m128i const & b) {
      __m128i mask_lo = _mm_load_si128((__m128i*)simd_mask_b);

      //== first shuffle (reverse) the bytes
      __m128i r = _mm_shuffle_epi8(b, _mm_load_si128((__m128i*)simd_rev_mask_b));

      // lower 4 bits shift to upper
      __m128i lo = _mm_slli_epi32(_mm_and_si128(mask_lo, r), 4);
      // upper 4 bits shift to lower.
      __m128i hi = _mm_srli_epi32(_mm_andnot_si128(mask_lo, r), 4);

      // swap bits in groups of 4 bits, so after this point, just or.

      return _mm_or_si128(hi, lo);
    }
    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value, __m128i>::type word_reverse_complement_simd(__m128i const & b) {
      __m128i mask_lo = _mm_load_si128((__m128i*)simd_mask_b);

      //== first shuffle (reverse) the bytes
      __m128i r = _mm_shuffle_epi8(b, _mm_load_si128((__m128i*)simd_rev_mask_b));

      //== get the lut indices.
      // lower 4 bits
      __m128i lo = _mm_and_si128(mask_lo, r);
      // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
      __m128i hi = _mm_srli_epi32(_mm_andnot_si128(mask_lo, r), 4);

      //== now use shuffle to look up the reversed bytes.  hi and lo are now indices. and lut2_lo/hi are lookup tables. remember that lo and hi need to be swapped.
      return _mm_or_si128(_mm_shuffle_epi8(_mm_load_si128((__m128i*)simd_lut1_hi_b), lo),
                          _mm_shuffle_epi8(_mm_load_si128((__m128i*)simd_lut1_lo_b), hi));

    }
//
//    void print128_num(__m128i var)
//    {
//        uint8_t *val = (uint8_t*) &var;
//        printf("Numerical: %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x \n",
//               val[0], val[1], val[2], val[3], val[4], val[5],
//               val[6], val[7], val[8], val[9], val[10], val[11], val[12], val[13],
//               val[14], val[15]);
//    }
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                  ::std::is_same<A, bliss::common::RNA>::value ||
                                  ::std::is_same<A, bliss::common::DNA16>::value, int>::type = 0>
    INLINE Kmer reverse_simd(Kmer const & src) const
    {
      //std::cout << "SIMD before = " << src << std::endl;

      Kmer result;  // empty kmer for output

      // do the rem.  When rem is not 0, then it's very slow.  (like 40 to 50 times slower.)  (because of maskmoveu)
      if (simd_rem > 0) {

        // not enough room, so need to copy.
        if (Kmer::nWords < simd_stride) {
          // allocate space
          WORD_TYPE tmp[simd_stride] alignas(16);

          // copy then load to register.
          memcpy((tmp + simd_stride - Kmer::nWords), src.getConstData(), Kmer::nWords * sizeof(WORD_TYPE));
          __m128i in = _mm_loadu_si128((__m128i*)tmp);

//          print(in , "in");

          // do reverse.
          __m128i mword = word_reverse_simd<A>(in);     // SSE2

          // copy back out
          _mm_storeu_si128((__m128i*)tmp, mword);

          // maskmoveu is very slow.  why?
          //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.getData()));               // SSE2
          memcpy(result.getData(), tmp, Kmer::nWords * sizeof(WORD_TYPE));                                             // SSE2

        } else {
          // has room, extra will be overwritten, so just use _mm_storeu.
          __m128i in = _mm_loadu_si128((__m128i*)(src.getConstData() + Kmer::nWords - simd_stride));

          __m128i mword = word_reverse_simd<A>(in);     // SSE2

          // store back out
          _mm_storeu_si128((__m128i*)(result.getData()), mword);                                                     // SSE2
        }

      }  // else no rem.


      // swap the word order.  also swap the packed chars in the word at the same time.
      if (simd_iters > 0) {
        auto in = src.getConstData();
        auto out = result.getData() + Kmer::nWords - simd_stride;

        for (int i = 0; i < simd_iters; ++i) {
          // load, reverse, then store.  assuming unaligned.
          _mm_storeu_si128((__m128i*)out, word_reverse_simd<A>(_mm_loadu_si128((__m128i*)in)));                 // SSE2
          in += simd_stride;
          out -= simd_stride;
        }
      }

//
      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.right_shift_bits(Kmer::nWords * sizeof(WORD_TYPE) * 8 - Kmer::nBits);

//      std::cout << "SIMD reversed, shifted = " << result << std::endl;
//      std::cout << "SIMD before = " << src << std::endl << std::endl;

      return result;
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                  ::std::is_same<A, bliss::common::RNA>::value ||
                                  ::std::is_same<A, bliss::common::DNA16>::value, int>::type = 0>
    INLINE Kmer reverse_complement_simd(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output

      // do the rem.  When rem is not 0, then it's very slow.  (like 40 to 50 times slower.)  (because of maskmoveu)
      if (simd_rem > 0) {

        // not enough room, so need to copy.
        if (Kmer::nWords < simd_stride) {
          // allocate space
          WORD_TYPE tmp[simd_stride] alignas(16);

          // copy then load to register.
          memcpy((tmp + simd_stride - Kmer::nWords), src.getConstData(), Kmer::nWords * sizeof(WORD_TYPE));
          __m128i in = _mm_loadu_si128((__m128i*)tmp);

//          print(in , "in");

          // do reverse.
          __m128i mword = word_reverse_complement_simd<A>(in);     // SSE2

          // copy back out
          _mm_storeu_si128((__m128i*)tmp, mword);

          // maskmoveu is very slow.  why?
          //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.getData()));               // SSE2
          memcpy(result.getData(), tmp, Kmer::nWords * sizeof(WORD_TYPE));                                             // SSE2

        } else {
          // has room, extra will be overwritten, so just use _mm_storeu.
          __m128i in = _mm_loadu_si128((__m128i*)(src.getConstData() + Kmer::nWords - simd_stride));

          __m128i mword = word_reverse_complement_simd<A>(in);     // SSE2

          // store back out
          _mm_storeu_si128((__m128i*)(result.getData()), mword);                                                     // SSE2
        }
      }  // else no rem.


      // swap the word order.  also swap the packed chars in the word at the same time.
      if (simd_iters > 0) {
        auto in = src.getConstData();
        auto out = result.getData() + Kmer::nWords - simd_stride;

        for (int i = 0; i < simd_iters; ++i) {
          // load, reverse, then store.  assuming unaligned.
          _mm_storeu_si128((__m128i*)out, word_reverse_complement_simd<A>(_mm_loadu_si128((__m128i*)in)));                // SSE2
          in += simd_stride;
          out -= simd_stride;
        }
      }



      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.right_shift_bits(Kmer::nWords * sizeof(WORD_TYPE) * 8 - Kmer::nBits);

      return result;
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, bliss::common::DNA>::value ||
                                    ::std::is_same<A, bliss::common::RNA>::value ||
                                    ::std::is_same<A, bliss::common::DNA16>::value), int>::type = 0>
    INLINE Kmer reverse_simd(Kmer const & src) const
    {
       return reverse_serial(src);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, bliss::common::DNA>::value ||
                                    ::std::is_same<A, bliss::common::RNA>::value ||
                                    ::std::is_same<A, bliss::common::DNA16>::value), int>::type = 0>
    INLINE Kmer reverse_complement_simd(Kmer const & src) const
    {
      return reverse_complement_serial(src);
    }


#endif

    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits. 
    /// 8 to 100 times faster when using DNA/RNA, and 8 to 50 times faster than sequentially reverse the bits when DNA16.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<(::std::is_same<A, bliss::common::DNA>::value ||
                                     ::std::is_same<A, bliss::common::RNA>::value ||
                                     ::std::is_same<A, bliss::common::DNA16>::value) &&
                                      (sizeof(W) == 8), W>::type word_reverse(W const & b) const {
      W v = b;
      W mask = (static_cast<W>(~0) >> 32);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
      v =  (v >> 32)         |  (v         << 32);  mask ^= (mask << 16); // swap 32 bit. can use BSWAP_64 here
      v = ((v >> 16) & mask) | ((v & mask) << 16);  mask ^= (mask << 8);  // swap 16 bit. can use BSWAP here.
      v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);  // sqap bytes. can use XCHG here
      v = ((v >>  4) & mask) | ((v & mask) <<  4);                        // swap nibbles.
      if (!::std::is_same<A, bliss::common::DNA16>::value) {                             // if not DNA16, then swap 2 bits
        mask ^= (mask << 2);
        v = ((v >>  2) & mask) | ((v & mask) <<  2);
      }

      return v;
    }
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<(::std::is_same<A, bliss::common::DNA>::value ||
                                     ::std::is_same<A, bliss::common::RNA>::value ||
                                     ::std::is_same<A, bliss::common::DNA16>::value) &&
                                      (sizeof(W) == 4), W>::type word_reverse(W const & b) const {
      W v = b;
      W mask = (static_cast<W>(~0) >> 16);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
      v =  (v >> 16)         |  (v         << 16);  mask ^= (mask << 8);  // swap 16 bit. can use BSWAP here.
      v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);  // sqap bytes. can use XCHG here
      v = ((v >>  4) & mask) | ((v & mask) <<  4);                        // swap nibbles.
      if (!::std::is_same<A, bliss::common::DNA16>::value) {                             // if not DNA16, then swap 2 bits
        mask ^= (mask << 2);
        v = ((v >>  2) & mask) | ((v & mask) <<  2);
      }

      return v;
    }
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<(::std::is_same<A, bliss::common::DNA>::value ||
                                     ::std::is_same<A, bliss::common::RNA>::value ||
                                     ::std::is_same<A, bliss::common::DNA16>::value) &&
                                      (sizeof(W) == 2), W>::type word_reverse(W const & b) const {
      W v = b;
      W mask = (static_cast<W>(~0) >> 8);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
      v =  (v >>  8)         |  (v         <<  8);  mask ^= (mask << 4);  // sqap bytes. can use XCHG here
      v = ((v >>  4) & mask) | ((v & mask) <<  4);                        // swap nibbles.
      if (!::std::is_same<A, bliss::common::DNA16>::value) {                             // if not DNA16, then swap 2 bits
        mask ^= (mask << 2);
        v = ((v >>  2) & mask) | ((v & mask) <<  2);
      }

      return v;
    }
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<(::std::is_same<A, bliss::common::DNA>::value ||
                                     ::std::is_same<A, bliss::common::RNA>::value ||
                                     ::std::is_same<A, bliss::common::DNA16>::value) &&
                                      (sizeof(W) == 1), W>::type word_reverse(W const & b) const {
      W v = b;
      W mask = (static_cast<W>(~0) >> 4);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
      v =  (v >>  4)         |  (v         <<  4);                        // swap nibbles.
      if (!::std::is_same<A, bliss::common::DNA16>::value) {                             // if not DNA16, then swap 2 bits
        mask ^= (mask << 2);
        v = ((v >>  2) & mask) | ((v & mask) <<  2);
      }

      return v;
    }
    
    /// do reverse complement.  8 to 100 times faster than serially reversing the bits.
    template <typename A = ALPHABET>
        INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                         ::std::is_same<A, bliss::common::RNA>::value, WORD_TYPE>::type word_reverse_complement(WORD_TYPE const & b) const {
      // DNA type, requires that alphabet be setup so that negation produces the complement.
      return ~(word_reverse<A, WORD_TYPE>(b));
    }


    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value &&
            (sizeof(W) == 8), W>::type word_reverse_complement(W const & b) const {
      // DNA type:
      W v = b;
      W mask = (static_cast<W>(~0) >> 32);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.

      // REQUIRES that alphabet where complement is bit reversal of the original character.
      // since for DNA16, complement is bit reversal in a nibble, we do 2 extra iterations.
      v =  (v >> 32)         |  (v         << 32);  mask ^= (mask << 16);
      v = ((v >> 16) & mask) | ((v & mask) << 16);  mask ^= (mask << 8);
      v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);
      v = ((v >>  4) & mask) | ((v & mask) <<  4);  mask ^= (mask << 2);  // reverse nibbles within a byte
      v = ((v >>  2) & mask) | ((v & mask) <<  2);  mask ^= (mask << 1);  // last 2 steps create the complement
      v = ((v >>  1) & mask) | ((v & mask) <<  1);

      // if we were to walk through the data byte by byte, and do complement and swap the low and high 4 bits,
      // the code would be as SLOW as do_reverse_complement.

      return v;
    }
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value &&
            (sizeof(W) == 4), W>::type word_reverse_complement(W const & b) const {
      // DNA type:
      W v = b;
      W mask = (static_cast<W>(~0) >> 16);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.

      // REQUIRES that alphabet where complement is bit reversal of the original character.
      // since for DNA16, complement is bit reversal in a nibble, we do 2 extra iterations.
      v =  (v >> 16)         |  (v         << 16);  mask ^= (mask << 8);
      v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);
      v = ((v >>  4) & mask) | ((v & mask) <<  4);  mask ^= (mask << 2);  // reverse nibbles within a byte
      v = ((v >>  2) & mask) | ((v & mask) <<  2);  mask ^= (mask << 1);  // last 2 steps create the complement
      v = ((v >>  1) & mask) | ((v & mask) <<  1);

      // if we were to walk through the data byte by byte, and do complement and swap the low and high 4 bits,
      // the code would be as SLOW as do_reverse_complement.

      return v;
    }
    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value &&
            (sizeof(W) == 2), W>::type word_reverse_complement(W const & b) const {
      // DNA type:
      W v = b;
      W mask = (static_cast<W>(~0) >> 8);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.

      // REQUIRES that alphabet where complement is bit reversal of the original character.
      // since for DNA16, complement is bit reversal in a nibble, we do 2 extra iterations.
      v =  (v >>  8)         |  (v         <<  8);  mask ^= (mask << 4);
      v = ((v >>  4) & mask) | ((v & mask) <<  4);  mask ^= (mask << 2);  // reverse nibbles within a byte
      v = ((v >>  2) & mask) | ((v & mask) <<  2);  mask ^= (mask << 1);  // last 2 steps create the complement
      v = ((v >>  1) & mask) | ((v & mask) <<  1);

      // if we were to walk through the data byte by byte, and do complement and swap the low and high 4 bits,
      // the code would be as SLOW as do_reverse_complement.

      return v;
    }
    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value &&
            (sizeof(W) == 1), W>::type word_reverse_complement(W const & b) const {
      // DNA type:
      W v = b;
      W mask = (static_cast<W>(~0) >> 4);

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.

      // REQUIRES that alphabet where complement is bit reversal of the original character.
      // since for DNA16, complement is bit reversal in a nibble, we do 2 extra iterations.
      v =  (v >>  4)         |  (v         <<  4);  mask ^= (mask << 2);  // reverse nibbles within a byte
      v = ((v >>  2) & mask) | ((v & mask) <<  2);  mask ^= (mask << 1);  // last 2 steps create the complement
      v = ((v >>  1) & mask) | ((v & mask) <<  1);

      // if we were to walk through the data byte by byte, and do complement and swap the low and high 4 bits,
      // the code would be as SLOW as do_reverse_complement.

      return v;
    }
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                  ::std::is_same<A, bliss::common::RNA>::value ||
                                  ::std::is_same<A, bliss::common::DNA16>::value, int>::type = 0>
    INLINE Kmer reverse_swar(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output


      // swap the word order.  also swap the packed chars in the word at the same time.
      for (int i = 0, j = Kmer::nWords - 1, max = Kmer::nWords; i < max; ++i, --j)
        result.getData()[j] = word_reverse<A, WORD_TYPE>(src.getConstData()[i]);

//      std::cout << "SWAR reversed = " << result << std::endl;

      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.right_shift_bits(Kmer::nWords * sizeof(WORD_TYPE) * 8 - Kmer::nBits);

//      std::cout << "SWAR reversed, shifted = " << result << std::endl;
//      std::cout << "SWAR before: " << src << std::endl;

      return result;
    }
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                  ::std::is_same<A, bliss::common::RNA>::value ||
                                  ::std::is_same<A, bliss::common::DNA16>::value, int>::type = 0>
    INLINE Kmer reverse_complement_swar(Kmer const & src) const
    {
      Kmer result;  // empty kmerfor output

      // swap the word order.  also do rev_comp for the chars in the word at the same time.
      for (int i = 0, j = Kmer::nWords - 1, max = Kmer::nWords; i < max; ++i, --j)
        result.getData()[j] = word_reverse_complement<A>(src.getConstData()[i]);

      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.right_shift_bits(Kmer::nWords * sizeof(WORD_TYPE) * 8 - Kmer::nBits);

      return result;
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, bliss::common::DNA>::value ||
                                    ::std::is_same<A, bliss::common::RNA>::value ||
                                    ::std::is_same<A, bliss::common::DNA16>::value), int>::type = 0>
    INLINE Kmer reverse_swar(Kmer const & src) const
    {
       return reverse_serial(src);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, bliss::common::DNA>::value ||
                                    ::std::is_same<A, bliss::common::RNA>::value ||
                                    ::std::is_same<A, bliss::common::DNA16>::value), int>::type = 0>
    INLINE Kmer reverse_complement_swar(Kmer const & src) const
    {
      return reverse_complement_serial(src);
    }



    // ================ experiement with bswap operations.

    static constexpr uint64_t mask0 = 0x0F0F0F0F0F0F0F0F;
    static constexpr uint64_t mask1 = 0x3333333333333333;
    static constexpr uint64_t mask2 = 0x5555555555555555;
//    static constexpr uint32_t mask_32[3] = { 0x0F0F0F0F, 0x33333333, 0x55555555 };
//    static constexpr uint16_t mask_16[3] = { 0x0F0F, 0x3333, 0x5555 };
//    static constexpr uint8_t  mask_8[3]  = { 0x0F, 0x33, 0x55 };

    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.  8 to 100 times faster.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                     ::std::is_same<A, bliss::common::RNA>::value, W>::type word_reverse_bswap(W const & b) const {
      W v = b;
      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      switch (sizeof(W) * 8)
      {
        // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
        case 64 :  v = __bswap_64(v); break; // bswap_64 instruction
        case 32 :  v = __bswap_32(v); break; // bswap or bswap_32 instruction
        case 16 :  v = __bswap_16(v); break; // XCHG instruction
        default :  break;
      }
      // finally, bit swap within the bytes.  only need to go down to chars of 2 bits.
      v = ((v >>  4) & mask0) | ((v & mask0) <<  4);  // SWAR style
      v = ((v >>  2) & mask1) | ((v & mask1) <<  2);

      return v;
    }
    /// do reverse complement.  8 to 100 times faster than serially reversing the bits.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
        INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                         ::std::is_same<A, bliss::common::RNA>::value, W>::type word_reverse_complement_bswap(W const & b) const {
      // DNA type, requires that alphabet be setup so that negation produces the complement.
      return ~(word_reverse_bswap<A, W>(b));
    }

    // reverse via bit swapping in 4 bit increment.  this complement table is then used for lookup.
    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.
    /// 8 to 50 times faster than sequentially reverse the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
        INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value, W>::type  word_reverse_bswap(W const & b) const {
      W v = b;

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      switch (sizeof(W) * 8)
      {
        // use bswap_xx from gnu c library.
        case 64 :  v = __bswap_64(v); break;   // bswap_64 instruction           // LATER.
        case 32 :  v = __bswap_32(v); break;   // bswap or bswap_32 instruction  // PENTIUM
        case 16 :  v = __bswap_16(v); break;   // XCHG instruction
        default :  break;
      }
      // finally, bit swap within the bytes.  only need to go down to chars of 4 bits.
      v = ((v >>  4) & mask0) | ((v & mask0) <<  4);  // SWAR style

      return v;
    }
    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
        INLINE typename ::std::enable_if<::std::is_same<A, bliss::common::DNA16>::value, W>::type  word_reverse_complement_bswap(W const & b) const {
      W v = b;

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      switch (sizeof(W) * 8)
      {
        // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
        // use bswap_xx from gnu c library.
        case 64 :  v = __bswap_64(v); break;   // bswap_64 instruction
        case 32 :  v = __bswap_32(v); break;   // bswap or bswap_32 instruction
        case 16 :  v = __bswap_16(v); break;   // XCHG instruction
        default :  break;
      }

      // finally, bit swap within the bytes.  first swap the 4 bit characters
      v = ((v >>  4) & mask0) | ((v & mask0) <<  4);  // SWAR style

      // now the complement.
      // since for DNA16, complement is bit reversal in a nibble, we do 2 extra iterations.
      v = ((v >>  2) & mask1) | ((v & mask1) <<  2);
      v = ((v >>  1) & mask2) | ((v & mask2) <<  1);
      // REQUIRES that alphabet where complement is bit reversal of the original character.

      // NOTE: if on AMD, VPPERM can do 128 bit bit reverse,

      return v;
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                  ::std::is_same<A, bliss::common::RNA>::value ||
                                  ::std::is_same<A, bliss::common::DNA16>::value, int>::type = 0>
    INLINE Kmer reverse_bswap(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output

      // unchunked reverse works.
            for (int i = 0, j = Kmer::nWords - 1, max = Kmer::nWords; i < max; ++i, --j)
              result.getData()[j] = word_reverse_bswap<A, WORD_TYPE>(src.getConstData()[i]);


/*      //== chunked revcomp does NOT always work, although slightly faster, especially for small words.
      if (iters > 0) {
    	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
              const MACH_WORD_TYPE *s = reinterpret_cast<const MACH_WORD_TYPE*>(src.getConstData());
              WORD_TYPE *d = result.getData() + Kmer::nWords - stride;

              // swap the word order.  also swap the packed chars in the word at the same time.
              // do the first part of source in multiple of MACH_WORD_TYPE.
              for (int i = 0; i < iters; ++i, d -= stride)
                *(reinterpret_cast<MACH_WORD_TYPE*>(d)) = word_reverse_bswap<A, MACH_WORD_TYPE>(s[i]);
            }
      //      printf("1st step:\t\t%s\n", result.toAlphabetString().c_str());

            // do the rem.
            if (rem > 1) {
          	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
              MACH_WORD_TYPE mword = word_reverse_bswap<A, MACH_WORD_TYPE>(*(reinterpret_cast<const MACH_WORD_TYPE*>(src.getConstData() + Kmer::nWords - stride)));
              memcpy(result.getData(), &mword, rem * sizeof(WORD_TYPE));
            } else if (rem == 1) {
              result.getData()[0] = word_reverse_bswap<A, WORD_TYPE>(src.getConstData()[Kmer::nWords - 1]);
            }  // else no rem.
      //      printf("2nd step:\t\t%s\n", result.toAlphabetString().c_str());
*/

            // shift if necessary      // ununsed bits will be set to 0 by shift
      result.right_shift_bits(Kmer::nWords * sizeof(WORD_TYPE) * 8 - Kmer::nBits);

      return result;
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, bliss::common::DNA>::value ||
                                  ::std::is_same<A, bliss::common::RNA>::value ||
                                  ::std::is_same<A, bliss::common::DNA16>::value, int>::type = 0>
    INLINE Kmer reverse_complement_bswap(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output


      // swap the word order.  also do rev_comp for the chars in the word at the same time.
      for (int i = 0, j = Kmer::nWords - 1, max = Kmer::nWords; i < max; ++i, --j)
        result.getData()[j] = word_reverse_complement_bswap<A, WORD_TYPE>(src.getConstData()[i]);

/*      //== chunked revcomp does NOT always work, although slightly faster, especially for small words.
      if (iters > 0) {
    	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
        const MACH_WORD_TYPE *s = reinterpret_cast<const MACH_WORD_TYPE*>(src.getConstData());
        WORD_TYPE *d = result.getData() + Kmer::nWords - stride;

        // swap the word order.  also swap the packed chars in the word at the same time.
        // do the first part of source in multiple of MACH_WORD_TYPE.
        for (int i = 0; i < iters; ++i, d -= stride)
          *(reinterpret_cast<MACH_WORD_TYPE*>(d)) = word_reverse_complement_bswap<A, MACH_WORD_TYPE>(s[i]);
      }
//      printf("1st step:\t\t%s\n", result.toAlphabetString().c_str());

      // do the rem.
      if (rem > 1) {
    	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
        MACH_WORD_TYPE mword = word_reverse_complement_bswap<A, MACH_WORD_TYPE>(*(reinterpret_cast<const MACH_WORD_TYPE*>(src.getConstData() + Kmer::nWords - stride)));
        memcpy(result.getData(), &mword, rem * sizeof(WORD_TYPE));
      } else if (rem == 1) {
        result.getData()[0] = word_reverse_complement_bswap<A, WORD_TYPE>(src.getConstData()[Kmer::nWords - 1]);
      }  // else no rem.
//      printf("2nd step:\t\t%s\n", result.toAlphabetString().c_str());
*/
      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.right_shift_bits(Kmer::nWords * sizeof(WORD_TYPE) * 8 - Kmer::nBits);

      return result;
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, bliss::common::DNA>::value ||
                                    ::std::is_same<A, bliss::common::RNA>::value ||
                                    ::std::is_same<A, bliss::common::DNA16>::value), int>::type = 0>
    INLINE Kmer reverse_bswap(Kmer const & src) const
    {
       return reverse_serial(src);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, bliss::common::DNA>::value ||
                                    ::std::is_same<A, bliss::common::RNA>::value ||
                                    ::std::is_same<A, bliss::common::DNA16>::value), int>::type = 0>
    INLINE Kmer reverse_complement_bswap(Kmer const & src) const
    {
      return reverse_complement_serial(src);
    }


  };

  
  
  template<typename Kmer>
  constexpr uint8_t KmerReverseHelper<Kmer>::simd_mask_b[16];
  template<typename Kmer>
  constexpr uint8_t KmerReverseHelper<Kmer>::simd_rev_mask_b[16];
  template<typename Kmer>
  constexpr uint8_t KmerReverseHelper<Kmer>::simd_lut1_lo_b[16];
  template<typename Kmer>
  constexpr uint8_t KmerReverseHelper<Kmer>::simd_lut1_hi_b[16];
  template<typename Kmer>
  constexpr uint8_t KmerReverseHelper<Kmer>::simd_lut2_lo_b[16];
  template<typename Kmer>
  constexpr uint8_t KmerReverseHelper<Kmer>::simd_lut2_hi_b[16];

    } // namespace test
  } // namespace common
} // namespace bliss



#endif // BLISS_COMMON_KMER_H
