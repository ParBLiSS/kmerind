/**
 * @file    kmer_reverse.cpp
 * @ingroup
 * @author  tpan
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */

#include <string>
#include <random>

#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "utils/timer.hpp"
#include "iterators/transform_iterator.hpp"

template <unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
void benchmark_reverse() {

  using KMER = bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;
  KMER kmer;
  KMER rev, revcomp;
  KMER rev2, revcomp2;
  KMER rev3, revcomp3;
  KMER rev4, revcomp4;
  int iterations = 10000001;
  std::vector<KMER> kmers;

  printf("testing kmer with k=%d, nbits = %d, word bits =%lu. stride %d, iters %d, remainder %d, simd_stride %d, simd_iters %d, simd_remainder %d\n", KMER_SIZE, KMER::nBits, sizeof(WORD_TYPE) * 8, KMER::stride, KMER::iters, KMER::rem, KMER::simd_stride, KMER::simd_iters, KMER::simd_rem);

  TIMER_INIT(test);

  TIMER_START(test);


  srand(0);
  for (int i = 0; i < KMER_SIZE; ++i) {
    kmer.nextFromChar(rand() % ALPHABET::SIZE);
  }
  kmers.reserve(iterations);
  for (int i = 0; i  < iterations; ++i) {
    kmer.nextFromChar(rand() % ALPHABET::SIZE);

    kmers.push_back(KMER(kmer));
  }

  TIMER_END(test, "init kmers", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    rev ^= kmers[i].reverse_serial();
  }
  TIMER_END(test, "seq rev", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    revcomp ^= kmers[i].reverse_complement_serial();
  }
  TIMER_END(test, "seq revcomp", iterations);



  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    rev2 ^= kmers[i].reverse_swar();
  }
  TIMER_END(test, "swar rev", iterations);


  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    revcomp2 ^= kmers[i].reverse_complement_swar();
  }
  TIMER_END(test, "swar revcomp", iterations);



  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    rev3 ^= kmers[i].reverse_bswap();
  }
  TIMER_END(test, "bswap rev", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    revcomp3 ^= kmers[i].reverse_complement_bswap();
  }
  TIMER_END(test, "bswap revcomp", iterations);

#if defined(__SSSE3__)
  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    rev4 ^= kmers[i].reverse_simd();
  }
  TIMER_END(test, "ssse3 rev", iterations);

  TIMER_START(test);
  for (int i = 0; i  < iterations; ++i) {
    revcomp4 ^= kmers[i].reverse_complement_simd();
  }
  TIMER_END(test, "ssse3 revcomp", iterations);
#endif


  TIMER_REPORT(test, 0);

}


template <unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
void test_reverse() {

  using KMER = bliss::common::Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>;
  KMER kmer;
  KMER rev, revcomp;
  KMER rev2, revcomp2;
  KMER rev3, revcomp3;
  KMER rev4, revcomp4;
  int iterations = 10000001;

  printf("testing kmer with k=%d, nbits = %d, word bits =%lu. stride %d, iters %d, remainder %d, simd_stride %d, simd_iters %d, simd_remainder %d\n", KMER_SIZE, KMER::nBits, sizeof(WORD_TYPE) * 8, KMER::stride, KMER::iters, KMER::rem, KMER::simd_stride, KMER::simd_iters, KMER::simd_rem);

  srand(0);
  for (int i = 0; i < KMER_SIZE; ++i) {
    kmer.nextFromChar(rand() % ALPHABET::SIZE);
  }

  for (int i = 0; i  < iterations; ++i) {
    rev ^= kmer.reverse_serial();
    revcomp ^= kmer.reverse_complement_serial();

    rev2 ^= kmer.reverse_swar();
    revcomp2 ^= kmer.reverse_complement_swar();

    if (rev2 != rev) {
      printf("Input kmer %s iter %d\n", kmer.toAlphabetString().c_str(), i);
      printf("ERROR: swar reverse different:\n\trev\t%s\n\tswar\t%s\n", rev.toAlphabetString().c_str(), rev2.toAlphabetString().c_str());
      break;
    }

    if (revcomp2 != revcomp) {
      printf("Input kmer %s iter %d\n", kmer.toAlphabetString().c_str(),i);
      printf("ERROR: swar rev comp different:\n\trev\t%s\n\tswar\t%s\n", revcomp.toAlphabetString().c_str(), revcomp2.toAlphabetString().c_str());
      break;
    }

    rev3 ^= kmer.reverse_bswap();
    revcomp3 ^= kmer.reverse_complement_bswap();


    if (rev != rev3) {
      printf("Input kmer %s iter %d\n", kmer.toAlphabetString().c_str(),i);
      printf("ERROR: bswap reverse different:\n\trev\t%s\n\trev3\t%s\n", rev.toAlphabetString().c_str(), rev3.toAlphabetString().c_str());
      break;
    }

    if (revcomp != revcomp3) {
      printf("Input kmer %s iter %d\n", kmer.toAlphabetString().c_str(),i);
      printf("ERROR: bswap rev comp different:\n\trev\t%s\n\trev3\t%s\n", revcomp.toAlphabetString().c_str(), revcomp3.toAlphabetString().c_str());
      break;
    }
#if defined(__SSSE3__)
    rev4 ^= kmer.reverse_simd();
    revcomp4 ^= kmer.reverse_complement_simd();


    if (rev != rev4) {
      printf("Input kmer %s iter %d\n", kmer.toAlphabetString().c_str(),i);
      printf("ERROR: ssse3 reverse different:\n\trev\t%s\n\trev4\t%s\n", rev.toAlphabetString().c_str(), rev4.toAlphabetString().c_str());
      break;
    }

    if (revcomp != revcomp4) {
      printf("Input kmer %s iter %d\n", kmer.toAlphabetString().c_str(),i);
      printf("ERROR: ssse3 rev comp different:\n\trev\t%s\n\trev4\t%s\n", revcomp.toAlphabetString().c_str(), revcomp4.toAlphabetString().c_str());
      break;
    }
#endif

    kmer.nextFromChar(rand() % ALPHABET::SIZE);
  }

}

int main(int argc, char** argv) {
	int choice = 3;

	if (argc > 1) {
		if (strncmp(argv[1], "--simd-test", 11) == 0) choice = 1;
		else if (strncmp(argv[1], "--test", 6) == 0) choice = 2;
		// else same.
	}

if (choice == 1) {
#if defined(__SSSE3__)
  unsigned char mask4_d[16] __attribute__((aligned(16))) = {0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F,
                                                            0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F};

  __m128i mask4 = _mm_loadu_si128((__m128i*)mask4_d);

  unsigned char arr[16] __attribute__((aligned(16))) = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                                                        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};

  __m128i in = _mm_loadu_si128((__m128i*)arr);

  unsigned char out[16] __attribute__((aligned(16)));

  _mm_storeu_si128((__m128i*)out, in);

  if (!std::equal(out, out+16, arr)) {
    printf("not same\n");
  }



  unsigned char mask_d[32] __attribute__((aligned(16))) = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
                                                           0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
                                                           0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                                           0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

  __m128i mask;


  for (int i = 0; i < 16; ++i) {
    mask = _mm_loadu_si128((__m128i*)(mask_d + i));

    std::fill(out, out+16, 0);
    _mm_maskmoveu_si128(in, mask, (char *)out);

    if (!std::equal(out, out+16-i, arr)) {
      printf("selected not same for iter %d\n", i);
    }
    if ((i > 0) && !std::equal(out+16-i, out+16, mask_d + 16)) {
      printf("cleared not same for iter %d\n", i);

    }

  }

  __m128i lo = _mm_and_si128(mask4, in);
  _mm_storeu_si128((__m128i*)out, lo);
  if (!std::equal(out, out+16, arr)) {
    printf("lo and not same\n");
  }

  __m128i hi = _mm_andnot_si128(mask4, in);
  _mm_storeu_si128((__m128i*)out, hi);
  if (!std::equal(out, out+16, mask_d+16)) {
    printf("hi andnot not same\n");
    for (int i = 0; i < 16; ++i) {
      printf("i %d src %d out %d expected %d\n", i, arr[i], out[i], mask_d[i+16]);
    }
  }

  hi = _mm_srli_epi32(_mm_andnot_si128(mask4, in), 4);
  _mm_storeu_si128((__m128i*)out, hi);
  if (!std::equal(out, out+16, mask_d+16)) {
    printf("hi shift right not same\n");
  }

  lo = _mm_slli_epi32(_mm_and_si128(mask4, in), 4);
  _mm_storeu_si128((__m128i*)out, lo);
  bool same = true;
  for (int i = 0; i < 16; ++i) {
    same &= out[i] == (arr[i] << 4);
  }
  if (!same) {
    printf("lo shift left not same\n");
  }


  __m128i hl = _mm_or_si128(hi, lo);
  _mm_storeu_si128((__m128i*)out, hl);
  same = true;
  for (int i = 0; i < 16; ++i) {
    same &= (out[i] == (arr[i] << 4));

  }
  if (!same) {
    printf("lo|hi not same\n");
  }




  unsigned char smask_d[32] __attribute__((aligned(16))) = {0x0f, 0x0e, 0x0d, 0x0c, 0x0b, 0x0a, 0x09, 0x08,
                                                            0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00,
                                                            0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
                                                            0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  __m128i smask, shuf;
  same = true;

  for (int i = 0; i < 16; ++i) {
    smask = _mm_loadu_si128((__m128i*)(smask_d + i));

    shuf = _mm_shuffle_epi8(in, smask);
    _mm_storeu_si128((__m128i*)out, shuf);

    for (int j = 15-i, k = 0; j >=0; --j, ++k) {
      same &= (out[k] == arr[j]);
    }
    if (!same) {
      printf("shuffled not same for iter %d\n", i);
    }
    if (i > 0) {
      for (int k = 15-i; k < 16; ++k) {
        same &= (out[k] == 0);
      }
    }
    if (!same) {
      printf("cleared not same for iter %d\n", i);
    }

  }

//  using KM = bliss::common::Kmer<32, bliss::common::DNA16, uint64_t>;
//  unsigned char rev_dna16[16] __attribute__((aligned(16))) = {0xf0, 0xe0, 0xd0, 0xc0, 0xb0, 0xa0, 0x90, 0x80,
//                                                               0x70, 0x60, 0x50, 0x40, 0x30, 0x20, 0x10, 0x00};
//
//  __m128i rev = KM::word_reverse_simd<bliss::common::DNA16>(in);
//  _mm_storeu_si128((__m128i*)out, rev);
//  if (!std::equal(out, out + 16, rev_dna16)) {
//    printf("rev16 not same.\n");
//    for (int i = 0; i < 16; ++i) {
//      printf("i %d src %x out %x expected %x\n", i, arr[i], out[i], rev_dna16[i]);
//    }
//  }
//
//  unsigned char revcomp_dna16[16] __attribute__((aligned(16))) = {0xf0, 0x70, 0xb0, 0x30, 0xd0, 0x50, 0x90, 0x10,
//                                                                   0xe0, 0x60, 0xa0, 0x20, 0xc0, 0x40, 0x80, 0x00};
//
//
//  __m128i revcomp = KM::word_reverse_complement_simd<bliss::common::DNA16>(in);
//  _mm_storeu_si128((__m128i*)out, revcomp);
//  if (!std::equal(out, out + 16, revcomp_dna16)) {
//    printf("revcomp16 not same.\n");
//    for (int i = 0; i < 16; ++i) {
//      printf("i %d src %x out %x expected %x\n", i, arr[i], out[i], revcomp_dna16[i]);
//    }
//  }
//
//
//  unsigned char rev_dna[16] __attribute__((aligned(16))) = {0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20,
//                                                             0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x00};
//
//  rev = KM::word_reverse_simd<bliss::common::DNA>(in);
//  _mm_storeu_si128((__m128i*)out, rev);
//  if (!std::equal(out, out + 16, rev_dna)) {
//    printf("rev4 not same.\n");
//    for (int i = 0; i < 16; ++i) {
//      printf("i %d src %x out %x expected %x\n", i, arr[i], out[i], rev_dna[i]);
//    }
//  }
//
//
//  revcomp = KM::word_reverse_complement_simd<bliss::common::DNA>(in);
//  _mm_storeu_si128((__m128i*)out, revcomp);
//  same = true;
//  for (int i = 0; i < 16; ++i) {
//    same &= (out[i] == static_cast<uint8_t>(~(rev_dna[i])));
//  }
//  if (!same) {
//    printf("revcomp4 not same\n");
//    for (int i = 0; i < 16; ++i) {
//      printf("i %d src %x out %x expected %x\n", i, arr[i], out[i], static_cast<uint8_t>(~(rev_dna[i])));
//    }
//  }
#else
  printf("SSSE3 is not defined and is needed for this test\n");
#endif
} else if (choice == 2) {



  test_reverse< 31, bliss::common::DNA, uint64_t>();  // 1 word, not full
  test_reverse< 32, bliss::common::DNA, uint64_t>();  // 1 word, full
  test_reverse< 33, bliss::common::DNA, uint64_t>();  // 2 words, not full
  test_reverse< 64, bliss::common::DNA, uint64_t>();  // 2 words, full
  test_reverse< 80, bliss::common::DNA, uint64_t>();  // 3 words, not full
  test_reverse< 96, bliss::common::DNA, uint64_t>();  // 3 words, full

  test_reverse< 15, bliss::common::DNA, uint32_t>();  // 1 word, not full
  test_reverse< 16, bliss::common::DNA, uint32_t>();  // 1 word, full
  test_reverse< 17, bliss::common::DNA, uint32_t>();  // 2 words, not full
  test_reverse< 32, bliss::common::DNA, uint32_t>();  // 2 words, full
  test_reverse< 40, bliss::common::DNA, uint32_t>();  // 3 words, not full
  test_reverse< 48, bliss::common::DNA, uint32_t>();  // 3 words, full

  test_reverse< 7, bliss::common::DNA, uint16_t>();  // 1 word, not full
  test_reverse< 8, bliss::common::DNA, uint16_t>();  // 1 word, full
  test_reverse< 9, bliss::common::DNA, uint16_t>();  // 2 words, not full
  test_reverse< 16, bliss::common::DNA, uint16_t>();  // 2 words, full
  test_reverse< 20, bliss::common::DNA, uint16_t>();  // 3 words, not full
  test_reverse< 24, bliss::common::DNA, uint16_t>();  // 3 words, full

  test_reverse< 3,  bliss::common::DNA, uint8_t>();  // 1 word, not full
  test_reverse< 4,  bliss::common::DNA, uint8_t>();  // 1 word, full
  test_reverse< 5,  bliss::common::DNA, uint8_t>();  // 2 words, not full
  test_reverse< 8,  bliss::common::DNA, uint8_t>();  // 2 words, full
  test_reverse< 10, bliss::common::DNA, uint8_t>();  // 3 words, not full
  test_reverse< 12, bliss::common::DNA, uint8_t>();  // 3 words, full

//  test_reverse< 21, bliss::common::DNA5, uint64_t>();  // 1 word, not full
//  test_reverse< 22, bliss::common::DNA5, uint64_t>();  // 2 word, not full
//  test_reverse< 42, bliss::common::DNA5, uint64_t>();  // 2 words, not full
//  test_reverse< 43, bliss::common::DNA5, uint64_t>();  // 3 words, not full
//  test_reverse< 64, bliss::common::DNA5, uint64_t>();  // 3 words, full
//
//  test_reverse< 2,  bliss::common::DNA5, uint8_t>();  // 1 word, not full
//  test_reverse< 3,  bliss::common::DNA5, uint8_t>();  // 2 word, not full
//  test_reverse< 5,  bliss::common::DNA5, uint8_t>();  // 2 words, not full
//  test_reverse< 6,  bliss::common::DNA5, uint8_t>();  // 3 words, not full
//  test_reverse< 8,  bliss::common::DNA5, uint8_t>();  // 3 words, full

  test_reverse< 15, bliss::common::DNA16, uint64_t>();  // 1 word, not full
  test_reverse< 16, bliss::common::DNA16, uint64_t>();  // 1 word, full
  test_reverse< 17, bliss::common::DNA16, uint64_t>();  // 2 words, not full
  test_reverse< 32, bliss::common::DNA16, uint64_t>();  // 2 words, full
  test_reverse< 40, bliss::common::DNA16, uint64_t>();  // 3 words, not full
  test_reverse< 48, bliss::common::DNA16, uint64_t>();  // 3 words, full

  test_reverse< 7, bliss::common::DNA16, uint32_t>();  // 1 word, not full
  test_reverse< 8, bliss::common::DNA16, uint32_t>();  // 1 word, full
  test_reverse< 9, bliss::common::DNA16, uint32_t>();  // 2 words, not full
  test_reverse< 16, bliss::common::DNA16, uint32_t>();  // 2 words, full
  test_reverse< 20, bliss::common::DNA16, uint32_t>();  // 3 words, not full
  test_reverse< 24, bliss::common::DNA16, uint32_t>();  // 3 words, full

  test_reverse< 3, bliss::common::DNA16, uint16_t>();  // 1 word, not full
  test_reverse< 4, bliss::common::DNA16, uint16_t>();  // 1 word, full
  test_reverse< 5, bliss::common::DNA16, uint16_t>();  // 2 words, not full
  test_reverse< 8, bliss::common::DNA16, uint16_t>();  // 2 words, full
  test_reverse< 10, bliss::common::DNA16, uint16_t>();  // 3 words, not full
  test_reverse< 12, bliss::common::DNA16, uint16_t>();  // 3 words, full


  test_reverse< 1,  bliss::common::DNA16, uint8_t >();  // 1 word, not full
  test_reverse< 2,  bliss::common::DNA16, uint8_t >();  // 1 word, full
  test_reverse< 3,  bliss::common::DNA16, uint8_t >();  // 2 words, not full
  test_reverse< 4,  bliss::common::DNA16, uint8_t >();  // 2 words, full
  test_reverse< 5,  bliss::common::DNA16, uint8_t >();  // 3 words, not full
  test_reverse< 6,  bliss::common::DNA16, uint8_t >();  // 3 words, full


} else {


  benchmark_reverse< 31, bliss::common::DNA, uint64_t>();  // 1 word, not full
  benchmark_reverse< 32, bliss::common::DNA, uint64_t>();  // 1 word, full
  benchmark_reverse< 33, bliss::common::DNA, uint64_t>();  // 2 words, not full
  benchmark_reverse< 64, bliss::common::DNA, uint64_t>();  // 2 words, full
  benchmark_reverse< 80, bliss::common::DNA, uint64_t>();  // 3 words, not full
  benchmark_reverse< 96, bliss::common::DNA, uint64_t>();  // 3 words, full

  benchmark_reverse< 31, bliss::common::DNA, uint32_t>();  // 1 word, not full
  benchmark_reverse< 32, bliss::common::DNA, uint32_t>();  // 1 word, full
  benchmark_reverse< 33, bliss::common::DNA, uint32_t>();  // 2 words, not full
  benchmark_reverse< 64, bliss::common::DNA, uint32_t>();  // 2 words, full
  benchmark_reverse< 80, bliss::common::DNA, uint32_t>();  // 3 words, not full
  benchmark_reverse< 96, bliss::common::DNA, uint32_t>();  // 3 words, full

  benchmark_reverse< 31, bliss::common::DNA, uint16_t>();  // 1 word, not full
  benchmark_reverse< 32, bliss::common::DNA, uint16_t>();  // 1 word, full
  benchmark_reverse< 33, bliss::common::DNA, uint16_t>();  // 2 words, not full
  benchmark_reverse< 64, bliss::common::DNA, uint16_t>();  // 2 words, full
  benchmark_reverse< 80, bliss::common::DNA, uint16_t>();  // 3 words, not full
  benchmark_reverse< 96, bliss::common::DNA, uint16_t>();  // 3 words, full

  benchmark_reverse< 3,  bliss::common::DNA, uint8_t>();  // 1 word, not full
  benchmark_reverse< 4,  bliss::common::DNA, uint8_t>();  // 1 word, full
  benchmark_reverse< 5,  bliss::common::DNA, uint8_t>();  // 2 words, not full
  benchmark_reverse< 8,  bliss::common::DNA, uint8_t>();  // 2 words, full
  benchmark_reverse< 10, bliss::common::DNA, uint8_t>();  // 3 words, not full
  benchmark_reverse< 12, bliss::common::DNA, uint8_t>();  // 3 words, full

//  benchmark_reverse< 21, bliss::common::DNA5, uint64_t>();  // 1 word, not full
//  benchmark_reverse< 22, bliss::common::DNA5, uint64_t>();  // 2 word, not full
//  benchmark_reverse< 42, bliss::common::DNA5, uint64_t>();  // 2 words, not full
//  benchmark_reverse< 43, bliss::common::DNA5, uint64_t>();  // 3 words, not full
//  benchmark_reverse< 64, bliss::common::DNA5, uint64_t>();  // 3 words, full
//
//  benchmark_reverse< 2,  bliss::common::DNA5, uint8_t>();  // 1 word, not full
//  benchmark_reverse< 3,  bliss::common::DNA5, uint8_t>();  // 2 word, not full
//  benchmark_reverse< 5,  bliss::common::DNA5, uint8_t>();  // 2 words, not full
//  benchmark_reverse< 6,  bliss::common::DNA5, uint8_t>();  // 3 words, not full
//  benchmark_reverse< 8,  bliss::common::DNA5, uint8_t>();  // 3 words, full

  benchmark_reverse< 15, bliss::common::DNA16, uint64_t>();  // 1 word, not full
  benchmark_reverse< 16, bliss::common::DNA16, uint64_t>();  // 1 word, full
  benchmark_reverse< 17, bliss::common::DNA16, uint64_t>();  // 2 words, not full
  benchmark_reverse< 32, bliss::common::DNA16, uint64_t>();  // 2 words, full
  benchmark_reverse< 40, bliss::common::DNA16, uint64_t>();  // 3 words, not full
  benchmark_reverse< 48, bliss::common::DNA16, uint64_t>();  // 3 words, full

  benchmark_reverse< 15, bliss::common::DNA16, uint32_t>();  // 1 word, not full
  benchmark_reverse< 16, bliss::common::DNA16, uint32_t>();  // 1 word, full
  benchmark_reverse< 17, bliss::common::DNA16, uint32_t>();  // 2 words, not full
  benchmark_reverse< 32, bliss::common::DNA16, uint32_t>();  // 2 words, full
  benchmark_reverse< 40, bliss::common::DNA16, uint32_t>();  // 3 words, not full
  benchmark_reverse< 48, bliss::common::DNA16, uint32_t>();  // 3 words, full

  benchmark_reverse< 15, bliss::common::DNA16, uint16_t>();  // 1 word, not full
  benchmark_reverse< 16, bliss::common::DNA16, uint16_t>();  // 1 word, full
  benchmark_reverse< 17, bliss::common::DNA16, uint16_t>();  // 2 words, not full
  benchmark_reverse< 32, bliss::common::DNA16, uint16_t>();  // 2 words, full
  benchmark_reverse< 40, bliss::common::DNA16, uint16_t>();  // 3 words, not full
  benchmark_reverse< 48, bliss::common::DNA16, uint16_t>();  // 3 words, full


  benchmark_reverse< 1,  bliss::common::DNA16, uint8_t >();  // 1 word, not full
  benchmark_reverse< 2,  bliss::common::DNA16, uint8_t >();  // 1 word, full
  benchmark_reverse< 3,  bliss::common::DNA16, uint8_t >();  // 2 words, not full
  benchmark_reverse< 4,  bliss::common::DNA16, uint8_t >();  // 2 words, full
  benchmark_reverse< 5,  bliss::common::DNA16, uint8_t >();  // 3 words, not full
  benchmark_reverse< 6,  bliss::common::DNA16, uint8_t >();  // 3 words, full
}
}
