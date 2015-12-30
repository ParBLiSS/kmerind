// include google test
#include <gtest/gtest.h>
//#include <boost/concept_check.hpp>

// include classes to test
#include "common/kmer.hpp"
#include "common/alphabets.hpp"
#include "utils/kmer_utils.hpp"
#include "iterators/transform_iterator.hpp"
#include "utils/logging.h"

template<unsigned int KMER_SIZE, typename ALPHABET, typename word_type=WordType>
using MyKmer = bliss::common::Kmer<KMER_SIZE, ALPHABET, word_type>;

// templated test function
template<typename kmer_word_type, typename input_word_type, unsigned int kmer_size=31, class Alphabet=bliss::common::DNA>
void test_kmer_with_word_type_packed(input_word_type* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned step=1) {

  typedef MyKmer<kmer_size, Alphabet, kmer_word_type> kmer_type;

  // create (fill) Kmer
  kmer_type kmer;

  // the expected value is a 64bit kmer.  we only need the prefix corresponding to the kmer_size.
  constexpr size_t expected_shift = ((64 / bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar() - kmer_size) * bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar());
  //BL_INFOF("expected_shift: %lu\n", expected_shift);

  input_word_type* kmer_pointer = kmer_data;
  //BL_INFOF("kmer pointer: %X\n", *kmer_pointer);

  // fill first kmer
  unsigned int offset = 0;
  offset = kmer.fillFromPackedStream(kmer_pointer, offset);

  uint64_t expected = *kmer_ex >> expected_shift;
  kmer_type kmer_ex_0(reinterpret_cast<kmer_word_type*>(&expected));
  //BL_INFOF("iter i = %d, expected: %016lX, actual %016lX\n", 0, kmer_ex_0.getPrefix64(), kmer.getPrefix64());

  EXPECT_EQ(kmer_ex_0, kmer) << "Kmer from stream should be equal to kmer from non-stream";


  // generate more kmers
  for (unsigned int i = step; i < nkmers; i += step)
  {
    kmer.nextFromPackedStream(kmer_pointer, offset);
    expected = kmer_ex[i] >> expected_shift;
    kmer_type kmer_ex_i(reinterpret_cast<kmer_word_type*>(&expected));

    //BL_INFOF("iter i = %d, expected: %016lX, actual %016lX\n", i, kmer_ex_i.getPrefix64(), kmer.getPrefix64());

    EXPECT_EQ(kmer_ex_i, kmer) << "Kmer compare unequal for sizeof(input)="<< sizeof(input_word_type) << ", sizeof(kmer_word)=" << sizeof(kmer_word_type) << ", size=" << kmer_size << ", bits=" << bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar() << " i = " << i;
  }
}



template<typename input_word_type, unsigned int kmer_size=31, class Alphabet=bliss::common::DNA>
void test_kmers_with_packed_input(input_word_type* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned int step=1)
{
  // test with various kmer base types
  test_kmer_with_word_type_packed<uint8_t,  input_word_type, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_packed<uint16_t, input_word_type, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_packed<uint32_t, input_word_type, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_packed<uint64_t, input_word_type, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
}

template<typename input_word_type>
void test_kmers_packed(input_word_type* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // test for bits per character: 2, 4, and 8 (no padding only!)
  test_kmers_with_packed_input<input_word_type, 31, bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 28, bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 13, bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 4,  bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 1,  bliss::common::DNA>(kmer_data, kmer_ex, nkmers);

//  test_kmers_with_packed_input<input_word_type, 10, Bits4>(kmer_data, kmer_ex, nkmers, 2);
//  test_kmers_with_packed_input<input_word_type, 13, Bits4>(kmer_data, kmer_ex, nkmers, 2);
//
//  test_kmers_with_packed_input<input_word_type, 7, Bits8>(kmer_data, kmer_ex, nkmers, 4);
//  test_kmers_with_packed_input<input_word_type, 5, Bits8>(kmer_data, kmer_ex, nkmers, 4);
//
}

template<typename input_word_type>
void test_kmers_3_packed(input_word_type* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // maximum in 64 bits is 21
  test_kmers_with_packed_input<input_word_type, 21, bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 20, bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 13, bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 9,  bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<input_word_type, 1,  bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerGeneration, TestKmerGenerationPacked2)
{
  // test sequence: 0xabbacafebabe1234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabbacafebabe1234deadbeef01c0ffee
   * print(",\n".join([" "*24 + "0x" + hex(val << 2*i)[-33 : -17] for i in range(0,33)]))
   */
    uint64_t kmer_ex[33] = {0xabbacafebabe1234, 0xaeeb2bfaeaf848d3,
                            0xbbacafebabe1234d, 0xeeb2bfaeaf848d37,
                            0xbacafebabe1234de, 0xeb2bfaeaf848d37a,
                            0xacafebabe1234dea, 0xb2bfaeaf848d37ab,
                            0xcafebabe1234dead, 0x2bfaeaf848d37ab6,
                            0xafebabe1234deadb, 0xbfaeaf848d37ab6f,
                            0xfebabe1234deadbe, 0xfaeaf848d37ab6fb,
                            0xebabe1234deadbee, 0xaeaf848d37ab6fbb,
                            0xbabe1234deadbeef, 0xeaf848d37ab6fbbc,
                            0xabe1234deadbeef0, 0xaf848d37ab6fbbc0,
                            0xbe1234deadbeef01, 0xf848d37ab6fbbc07,
                            0xe1234deadbeef01c, 0x848d37ab6fbbc070,
                            0x1234deadbeef01c0, 0x48d37ab6fbbc0703,
                            0x234deadbeef01c0f, 0x8d37ab6fbbc0703f,
                            0x34deadbeef01c0ff, 0xd37ab6fbbc0703ff,
                            0x4deadbeef01c0ffe, 0x37ab6fbbc0703ffb,
                            0xdeadbeef01c0ffee};

  // unpadded stream (bits_per_char is 2 => no padding)
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*2)) & 0x3) << (j*2) ) for j in range(0,4)], 0) )[:-1] for i in range(0, 64, 4) ) )
   */
  // 8 bit input
  uint8_t kmer_data_8[16] = {
                           0xea,
                           0xae,
                           0xa3,
                           0xbf,
                           0xae,
                           0xbe,
                           0x84,
                           0x1c,
                           0xb7,
                           0x7a,
                           0xbe,
                           0xfb,
                           0x40,
                           0x03,
                           0xff,
                           0xbb
  };
  test_kmers_packed<uint8_t>(kmer_data_8, kmer_ex, 33);

  // 16 bit input.   128 bits in input.  this covers the initial 64 bits, and 32 2-bit shifts, for total of 33 entries.
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*2)) & 0x3) << (j*2) ) for j in range(0,8)], 0) )[:-1] for i in range(0, 64, 8) ) )
    */
  uint16_t kmer_data_16[8] = {
                              0xaeea,
                              0xbfa3,
                              0xbeae,
                              0x1c84,
                              0x7ab7,
                              0xfbbe,
                              0x340,
                              0xbbff
  };
  test_kmers_packed<uint16_t>(kmer_data_16, kmer_ex, 33);


  // 32 bit input.   128 bits in input.  this covers the initial 64 bits, and 32 2-bit shifts, for total of 33 entries.
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*2)) & 0x3) << (j*2) ) for j in range(0,16)], 0) )[:-1] for i in range(0, 64, 16) ) )
    */
  uint32_t kmer_data_32[4] = {
                              0xbfa3aeea,
                              0x1c84beae,
                              0xfbbe7ab7,
                              0xbbff0340
  };
  test_kmers_packed<uint32_t>(kmer_data_32, kmer_ex, 33);


  // 64 bit input.   128 bits in input.  this covers the initial 64 bits, and 32 2-bit shifts, for total of 33 entries.
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*2)) & 0x3) << (j*2) ) for j in range(0,32)], 0) )[:-1] for i in range(0, 64, 32) ) )
    */
  uint64_t kmer_data_64[2] = {
                              0x1c84beaebfa3aeea,
                              0xbbff0340fbbe7ab7
  };
  test_kmers_packed<uint64_t>(kmer_data_64, kmer_ex, 33);


}

/**
 * Test k-mer generation with 3 bits and thus padded input
 */
TEST(KmerGeneration, TestKmerGenerationPacked3)
{
  // test sequence: 0xabbacafebabe1234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabbacafebabe1234deadbeef01c0ffee
   * print(",\n".join([" "*24 + hex(val >> (128 - 63 - 3*i) & 0x7fffffffffffffff)[:-1] for i in range(0,22)]))
   */
  uint64_t kmer_ex[22] = {
                           0x55dd657f5d5f091a,
                           0x2eeb2bfaeaf848d3,
                           0x77595fd757c2469b,
                           0x3acafebabe1234de,
                           0x5657f5d5f091a6f5,
                           0x32bfaeaf848d37ab,
                           0x15fd757c2469bd5b,
                           0x2febabe1234deadb,
                           0x7f5d5f091a6f56df,
                           0x7aeaf848d37ab6fb,
                           0x5757c2469bd5b7dd,
                           0x3abe1234deadbeef,
                           0x55f091a6f56df778,
                           0x2f848d37ab6fbbc0,
                           0x7c2469bd5b7dde03,
                           0x61234deadbeef01c,
                           0x091a6f56df7780e0,
                           0x48d37ab6fbbc0703,
                           0x469bd5b7dde0381f,
                           0x34deadbeef01c0ff,
                           0x26f56df7780e07ff,
                           0x37ab6fbbc0703ffb
  };
  // 8 bit input
  /* python:
    *  ... print(",\n".join(hex( (((val >> (128 - (i+2) * 3)) & 0x7) << 3) | ((val >> (128 - (i+1) * 3)) & 0x7) )[:-1] for i in range(0, 42, 2)))
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*3)) & 0x7) << (j*3) ) for j in range(0,2)], 0) )[:-1] for i in range(0, 42, 2) ) )

    */
 uint8_t kmer_data_8[21] = {
                          0x15,
                          0x1f,
                          0x1d,
                          0x11,
                          0x3f,
                          0x1d,
                          0x15,
                          0x37,
                          0x20,
                          0x1c,
                          0x1a,
                          0x33,
                          0x1d,
                          0x1b,
                          0x1f,
                          0x3d,
                          0x0,
                          0x23,
                          0x18,
                          0x3f,
                          0x1f
 };

  // unpadded stream (bits_per_char is 3, 2 chars per char => 2 bit padding)
  test_kmers_3_packed<uint8_t>(kmer_data_8, kmer_ex, 22);




  // 16 bit input  1 bit pad.   8 entries means only have 120 bits in input.  this covers the initial 63 bits, and 19 3-bit shifts, for total of 20 entries.
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*3)) & 0x7) << (j*3) ) for j in range(0,5)], 0) )[:-1] for i in range(0, 40, 5) ) )
    */
  uint16_t kmer_data_16[8] = {
                              0x57d5,
                              0x7e8b,
                              0x755d,
                              0x3906,
                              0x5cda,
                              0x3edb,
                              0x303d,
                              0x7ec4
  };

  // unpadded stream (bits_per_char is 3, 2 chars per char => 2 bit padding)
  test_kmers_3_packed<uint16_t>(kmer_data_16, kmer_ex, 20);


  // 32 bit input  2 bit pad.   4 entries means only have 120 bits in input.  this covers the initial 63 bits, and 19 3-bit shifts, for total of 20 entries.
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*3)) & 0x7) << (j*3) ) for j in range(0,10)], 0) )[:-1] for i in range(0, 40, 10) ) )
    */
  uint32_t kmer_data_32[4] = {
                              0x3f45d7d5,
                              0x1c83755d,
                              0x1f6ddcda,
                              0x3f62303d
  };

  // unpadded stream (bits_per_char is 3, 2 chars per char => 2 bit padding)
  test_kmers_3_packed<uint32_t>(kmer_data_32, kmer_ex, 20);

  // 64 bit input  1 bit pad.   2 entries means only have 126 bits in input.  this covers the initial 63 bits, and 21 3-bit shifts, for total of 22 entries.
  /* python:
   * import operator
   * print(",\n".join( hex( reduce(operator.or_, [( ((val >> (128 - (i+j+1)*3)) & 0x7) << (j*3) ) for j in range(0,21)], 0) )[:-1] for i in range(0, 42, 21) ) )
    */
  uint64_t kmer_data_64[2] = {
                              0x2720dd577f45d7d5,
                              0x3ffb1181ebedbb9b
  };

  // unpadded stream (bits_per_char is 3, 2 chars per char => 2 bit padding)
  test_kmers_3_packed<uint64_t>(kmer_data_64, kmer_ex, 22);


}



// templated test function
template<typename kmer_word_type, unsigned int kmer_size=31, class Alphabet=bliss::common::DNA>
void test_kmer_with_word_type_unpacked(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned step=1) {

  typedef MyKmer<kmer_size, Alphabet, kmer_word_type> kmer_type;

  // create (fill) Kmer
  kmer_type kmer;

  // the expected value is a 64bit kmer.  we only need the prefix corresponding to the kmer_size.
  constexpr size_t expected_shift = ((64 / bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar() - kmer_size) * bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar());
  //BL_INFOF("expected_shift: %lu\n", expected_shift);

  unsigned char* kmer_pointer = kmer_data;
  //BL_INFOF("kmer pointer: %d\n", *kmer_pointer);
  // fill first kmer
  //unsigned int offset = 0;
  kmer.fillFromChars(kmer_pointer);
  uint64_t expected = *kmer_ex >> expected_shift;
  kmer_type kmer_ex_0(reinterpret_cast<kmer_word_type*>(&expected));
  //BL_INFOF("iter i = %d, expected: %016lX, actual %016lX\n", 0, kmer_ex_0.getPrefix64(), kmer.getPrefix64());

  EXPECT_EQ(kmer, kmer_ex_0) << "Kmer from stream should be equal to kmer from non-stream";


  // generate more kmers
  for (unsigned int i = step; i < nkmers; i += step)
  {
    kmer.nextFromChar(*kmer_pointer); ++kmer_pointer;
    expected = kmer_ex[i] >> expected_shift;
    kmer_type kmer_ex_i(reinterpret_cast<kmer_word_type*>(&expected));

    //BL_INFOF("iter i = %d, expected: %016lX, actual %016lX\n", i, kmer_ex_i.getPrefix64(), kmer.getPrefix64());

    EXPECT_EQ(kmer_ex_i, kmer) << "Kmer compare unequal for sizeof(input)="<< sizeof(unsigned char) << ", sizeof(kmer_word)=" << sizeof(kmer_word_type) << ", size=" << kmer_size << ", bits=" << bliss::common::AlphabetTraits<Alphabet>::getBitsPerChar() << " i = " << i;
  }
}




template<unsigned int kmer_size=31, class Alphabet=bliss::common::DNA>
void test_kmers_with_unpacked_input(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned int step=1)
{
  // test with various kmer base types
  test_kmer_with_word_type_unpacked<uint8_t,  kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_unpacked<uint16_t, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_unpacked<uint32_t, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_unpacked<uint64_t, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
}

void test_kmers_unpacked(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // test for bits per character: 2, 4, and 8 (no padding only!)
  test_kmers_with_unpacked_input<31, bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<28, bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<13, bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<4,  bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<1,  bliss::common::DNA>(kmer_data, kmer_ex, nkmers);
}

void test_kmers_3_unpacked(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // maximum in 64 bits is 21
  test_kmers_with_unpacked_input<21, bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<20, bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<13, bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<9,  bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_unpacked_input<1,  bliss::common::DNA5>(kmer_data, kmer_ex, nkmers);
}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerGeneration, TestKmerGenerationChar2)
{
  // test sequence: 0xabbacafebabe1234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabbacafebabe1234deadbeef01c0ffee
   * print(",\n".join([" "*24 + "0x" + hex(val << 2*i)[-33 : -17] for i in range(0,33)]))
   */

    uint64_t kmer_ex[33] = {0xabbacafebabe1234, 0xaeeb2bfaeaf848d3,
                            0xbbacafebabe1234d, 0xeeb2bfaeaf848d37,
                            0xbacafebabe1234de, 0xeb2bfaeaf848d37a,
                            0xacafebabe1234dea, 0xb2bfaeaf848d37ab,
                            0xcafebabe1234dead, 0x2bfaeaf848d37ab6,
                            0xafebabe1234deadb, 0xbfaeaf848d37ab6f,
                            0xfebabe1234deadbe, 0xfaeaf848d37ab6fb,
                            0xebabe1234deadbee, 0xaeaf848d37ab6fbb,
                            0xbabe1234deadbeef, 0xeaf848d37ab6fbbc,
                            0xabe1234deadbeef0, 0xaf848d37ab6fbbc0,
                            0xbe1234deadbeef01, 0xf848d37ab6fbbc07,
                            0xe1234deadbeef01c, 0x848d37ab6fbbc070,
                            0x1234deadbeef01c0, 0x48d37ab6fbbc0703,
                            0x234deadbeef01c0f, 0x8d37ab6fbbc0703f,
                            0x34deadbeef01c0ff, 0xd37ab6fbbc0703ff,
                            0x4deadbeef01c0ffe, 0x37ab6fbbc0703ffb,
                            0xdeadbeef01c0ffee};


  // unpadded stream (bits_per_char is 2 => no padding)
  /* python:
   * print(",\n".join(hex((val >> (126 - i* 2)) & 0x3)[:-1] for i in range(0, 64)))
   */
  unsigned char kmer_data[64] = {2, 2, 2, 3, 2, 3, 2, 2,
                                 3, 0, 2, 2, 3, 3, 3, 2,
                                 2, 3, 2, 2, 2, 3, 3, 2,
                                 0, 1, 0, 2, 0, 3, 1, 0,
                                 3, 1, 3, 2, 2, 2, 3, 1,
                                 2, 3, 3, 2, 3, 2, 3, 3,
                                 0, 0, 0, 1, 3, 0, 0, 0,
                                 3, 3, 3, 3, 3, 2, 3, 2};

  // test with this data
  test_kmers_unpacked(kmer_data, kmer_ex, 33);

}


template<typename Alphabet, int K>
void compute_kmer(std::string input) {

  using KmerType = bliss::common::Kmer<K, Alphabet>;

  using Decoder = bliss::common::ASCII2<Alphabet, std::string::value_type>;
  using BaseCharIterator = bliss::iterator::transform_iterator<std::string::const_iterator, Decoder>;
  auto temp = BaseCharIterator(input.cbegin(), Decoder());

  KmerType kmer;
  kmer.fillFromChars(temp, false);
  unsigned int i = K;
  for (; i < input.length(); ++i) {
    std::string gold = input.substr(i-K, K);

    int res = strncmp(gold.c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str(), K);

    if (res != 0) {
      BL_INFOF("%d iterator input %s\n", i, gold.c_str());
      BL_INFOF("kmer %s %s %s\n", kmer.toString().c_str(), kmer.toAlphabetString().c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str());
    }

    EXPECT_EQ(res, 0);

    kmer.nextFromChar(*temp);  ++temp;

  }
  std::string gold = input.substr(input.length()-K, K);

  int res = strncmp(gold.c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str(), K);

  if (res != 0) {
    BL_INFOF("%d iterator input %s\n", i, gold.c_str());
    BL_INFOF("kmer %s %s %s\n", kmer.toString().c_str(), kmer.toAlphabetString().c_str(), bliss::utils::KmerUtils::toASCIIString(kmer).c_str());
  }

  EXPECT_EQ(res, 0);

}


/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerGeneration, TestKmerGenerationRoundTrip)
{
  // test sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";

  compute_kmer<bliss::common::DNA, 21>(input);
}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerGeneration, TestKmerGenerationRoundTripDNA5)
{
  // test sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";

  compute_kmer<bliss::common::DNA5, 21>(input);

}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerGeneration, TestKmerGenerationRoundTripMultiWord)
{
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";

  compute_kmer<bliss::common::DNA, 33>(input);

}

/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(KmerGeneration, TestKmerGenerationRoundTripDNA5MultiWord)
{
  std::string input = "GATTTGGGGTTCAAAGCAGT"
                         "ATCGATCAAATAGTAAATCC"
                         "ATTTGTTCAACTCACAGTTT";

  compute_kmer<bliss::common::DNA5, 33>(input);
}


/**
 * Test k-mer generation with 3 bits and thus padded input
 */
TEST(KmerGeneration, TestKmerGenerationChar3)
{
  // test sequence: 0xabbacafebabe1234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabbacafebabe1234deadbeef01c0ffee
   * print(",\n".join([" "*24 + hex(val >> (128 - 63 - 3*i) & 0x7fffffffffffffff)[:-1] for i in range(0,22)]))
   */
  uint64_t kmer_ex[22] = {
                           0x55dd657f5d5f091a,
                           0x2eeb2bfaeaf848d3,
                           0x77595fd757c2469b,
                           0x3acafebabe1234de,
                           0x5657f5d5f091a6f5,
                           0x32bfaeaf848d37ab,
                           0x15fd757c2469bd5b,
                           0x2febabe1234deadb,
                           0x7f5d5f091a6f56df,
                           0x7aeaf848d37ab6fb,
                           0x5757c2469bd5b7dd,
                           0x3abe1234deadbeef,
                           0x55f091a6f56df778,
                           0x2f848d37ab6fbbc0,
                           0x7c2469bd5b7dde03,
                           0x61234deadbeef01c,
                           0x091a6f56df7780e0,
                           0x48d37ab6fbbc0703,
                           0x469bd5b7dde0381f,
                           0x34deadbeef01c0ff,
                           0x26f56df7780e07ff,
                           0x37ab6fbbc0703ffb
  };

  /* python:
   *  print(",\n".join([" "*24 + hex(val >> (128 - 3 - 3*i) & 0x7)[:-1] for i in range(0,42)]))
   */
  uint8_t kmer_data_8[42] = {
                             0x5,
                             0x2,
                             0x7,
                             0x3,
                             0x5,
                             0x3,
                             0x1,
                             0x2,
                             0x7,
                             0x7,
                             0x5,
                             0x3,
                             0x5,
                             0x2,
                             0x7,
                             0x6,
                             0x0,
                             0x4,
                             0x4,
                             0x3,
                             0x2,
                             0x3,
                             0x3,
                             0x6,
                             0x5,
                             0x3,
                             0x3,
                             0x3,
                             0x7,
                             0x3,
                             0x5,
                             0x7,
                             0x0,
                             0x0,
                             0x3,
                             0x4,
                             0x0,
                             0x3,
                             0x7,
                             0x7,
                             0x7,
                             0x3
  };

  // test with 8 bit (padded by 2 bits)
  test_kmers_3_unpacked(kmer_data_8, kmer_ex, 22);

}


/**
 * Testing k-mer comparison operators.
 */
TEST(KmerComparison, TestKmerComparison1)
{
  // the main kmer value
  uint16_t kmer_val[] = {0xffee, 0x1c0, 0xbeef, 0xdead, 0x1234, 0x5678, 0xabba};
  // smaller value in 4th block::
  uint16_t kmer_val_s4[] = {0xffee, 0x1c0, 0xbeef, 0x1111, 0x1234, 0x5678, 0xabba};
  // greater value in 3rd block:
  uint16_t kmer_val_g3[] = {0xffee, 0x1c0, 0xfeef, 0xdead, 0x1234, 0x5678, 0xabba};

  MyKmer<41, bliss::common::DNA, uint16_t> kmer(kmer_val);
  MyKmer<41, bliss::common::DNA, uint16_t> kmer_s(kmer_val_s4);
  MyKmer<41, bliss::common::DNA, uint16_t> kmer_g(kmer_val_g3);

  EXPECT_TRUE(kmer > kmer_s);
  EXPECT_TRUE(kmer == kmer);
  EXPECT_TRUE(kmer_g > kmer);
  EXPECT_FALSE(kmer_g <= kmer);
  EXPECT_TRUE(kmer <= kmer);
  EXPECT_TRUE(kmer >= kmer);
  EXPECT_FALSE(kmer < kmer);
  EXPECT_FALSE(kmer > kmer);
  EXPECT_TRUE(kmer != kmer_g);
  EXPECT_TRUE(kmer != kmer_s);

}



/**
 * Testing kmer reverse
 */
TEST(KmerReverse, TestKmerReverse112)
{
  /*
   * python code to generate the reverse (by n bits each) for a given value:
   *
   * n = 2 # or n = 3
   * b = bin(val)
   * hex(int("".join(reversed([b[i+2+((len(b)-2) % n):][:n] for i in range(0, len(b)-2, n)])),2))
   */

  // testing with 112 bit sequence
  // test sequence: val = 0xabba56781234deadbeef01c0ffee
  // n = 2:
  // reverse seq:   val = 0xbbff0340fbbe7ab71c842d95aeea
  // n = 3
  // reverse seq:   val = 0x6bff23113ebedabd34a427952faa
  // n = 4
  // reverse seq:   val = 0xeeff0c10feebdaed43218765abba
  // n = 5
  // reverse seq:   val = 0x1dff8780e77cd5f5ba40b13ad375
  // n = 7
  // reverse seq:   val = 0xddfc18ee1777d6bda6440cf2b755

  /* python:
   * ", ".join(hex((val >> i*16) & 0xffff) for i in range(0, 8))
   */
  uint16_t kmer_val[] =  {0xffee, 0x1c0, 0xbeef, 0xdead, 0x1234, 0x5678, 0xabba};
  uint16_t kmer_ex[] =   {0xaeea, 0x2d95, 0x1c84, 0x7ab7, 0xfbbe, 0x340, 0xbbff};
  uint16_t kmer_ex_3[] = {0x2faa, 0x2795, 0x34a4, 0xdabd, 0x3ebe, 0x2311, 0x6bff};
  uint16_t kmer_ex_4[] = {0xabba, 0x8765, 0x4321, 0xdaed, 0xfeeb, 0xc10, 0xeeff};


  /* test for bits_per_char = 2 */
  MyKmer<7*8, bliss::common::DNA, uint16_t> kmer_in(kmer_val);
  MyKmer<7*8, bliss::common::DNA, uint16_t> kmer_ex_rev(kmer_ex);

//  std::cout << "src kmer: " << kmer_in << std::endl;

  MyKmer<7*8, bliss::common::DNA, uint16_t> kmer_rev = kmer_in.reverse();
//  std::cout << "rev kmer: " << kmer_rev << std::endl;
  EXPECT_EQ(kmer_ex_rev, kmer_rev);


  /* test for bits_per_char = 3 */
  MyKmer<37, bliss::common::DNA5, uint16_t> kmer3_in(kmer_val);
  MyKmer<37, bliss::common::DNA5, uint16_t> kmer3_ex_rev(kmer_ex_3);
  // get the reverse
  MyKmer<37, bliss::common::DNA5, uint16_t> kmer3_rev = kmer3_in.reverse();
  EXPECT_EQ(kmer3_ex_rev, kmer3_rev);

  /* test for bits_per_char = 4 */
  MyKmer<28, bliss::common::DNA16, uint16_t> kmer4_in(kmer_val);
  MyKmer<28, bliss::common::DNA16, uint16_t> kmer4_ex_rev(kmer_ex_4);
  // get the reverse
  MyKmer<28, bliss::common::DNA16, uint16_t> kmer4_rev = kmer4_in.reverse();
  EXPECT_EQ(kmer4_ex_rev, kmer4_rev);

}
