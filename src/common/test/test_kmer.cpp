// include google test
#include <gtest/gtest.h>
//#include <boost/concept_check.hpp>

// include classes to test
#include <common/Kmer.hpp>
#include <common/alphabets.hpp>


template<unsigned int KMER_SIZE, typename ALPHABET, typename word_type=WordType>
class MyKmer : bliss::Kmer<KMER_SIZE, ALPHABET, word_type>
{
  public:
    using SuperType = typename bliss::Kmer<KMER_SIZE, ALPHABET, word_type>;

    explicit MyKmer(const SuperType& other) : SuperType(other) {};

    template<typename InputIterator>
    explicit MyKmer(InputIterator iter) : SuperType(iter) {};
    MyKmer() : SuperType() {};
    template <typename InputIterator, typename offset_t>
    // TODO: add option for bit offset in input sequence?
    unsigned int fillFromPackedStream(InputIterator& begin, offset_t& offset, bool stop_on_last = false) {
      return SuperType::fillFromPackedStream(begin, offset, stop_on_last);
    }
    template <typename InputIterator>
    void fillFromChars(InputIterator& begin, bool stop_on_last = false) {
      SuperType::fillFromChars(begin, stop_on_last);
    }
    template <typename InputIterator, typename offset_t>
    void nextFromPackedStream(InputIterator& begin, offset_t& offset) {
      SuperType::nextFromPackedStream(begin, offset);
    }
    void nextFromChar(char c) {
      SuperType::nextFromChar(c);
    }

    inline bool operator==(const MyKmer& rhs) const
    {
      return SuperType::operator==(rhs);
    }
    inline bool operator!=(const MyKmer& rhs) const
    {
      return SuperType::operator!=(rhs);
    }
    inline bool operator<(const MyKmer& rhs) const
    {
      return SuperType::operator<(rhs);
    }
    inline bool operator<=(const MyKmer& rhs) const
    {
      return SuperType::operator<=(rhs);
    }
    inline bool operator>(const MyKmer& rhs) const
    {
      return SuperType::operator>(rhs);
    }
    inline bool operator>=(const MyKmer& rhs) const
    {
      return SuperType::operator>=(rhs);
    }
    inline MyKmer operator^(const MyKmer& rhs) const
    {
      return MyKmer(SuperType::operator^(rhs).data);
    }
    inline MyKmer& operator^=(const MyKmer& rhs) const
    {
      return SuperType::operator^=(rhs);
    }
    inline MyKmer operator&(const MyKmer& rhs) const
    {
      return SuperType::operator&(rhs);
    }
    inline MyKmer& operator&=(const MyKmer& rhs) const
    {
      return SuperType::operator&=(rhs);
    }
    inline MyKmer operator|(const MyKmer& rhs) const
    {
      return SuperType::operator|(rhs);
    }
    inline MyKmer& operator|=(const MyKmer& rhs) const
    {
      return SuperType::operator|=(rhs);
    }
    inline MyKmer operator<<(const MyKmer& rhs) const
    {
      return SuperType::operator<<(rhs);
    }
    inline MyKmer& operator<<=(const MyKmer& rhs) const
    {
      return SuperType::operator<<=(rhs);
    }
    inline MyKmer operator>>(const MyKmer& rhs) const
    {
      return SuperType::operator>>(rhs);
    }
    inline MyKmer& operator>>=(const MyKmer& rhs) const
    {
      return SuperType::operator>>=(rhs);
    }
    MyKmer reversed_kmer() const
    {
      return MyKmer(SuperType::reversed_kmer());
    }
    std::string toString() const
    {
      return SuperType::toString();
    }

    uint64_t getPrefix64() const {
      return SuperType::getPrefix64();
    }

};


struct Bits4
{
  static constexpr AlphabetSizeType SIZE = 15;
  // lookup table for XYZ
  static constexpr uint8_t FROM_ASCII[256] =
  {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'A'     'C'             'G'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 'T'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'a'     'c'             'g'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 't'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  };

  // reverse lookup table for XYZ
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T'  // = 3
  };

  // reverse lookup table for DNA5
  static constexpr uint8_t TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0  // = 3
  };
};
struct Bits5
{
  static constexpr AlphabetSizeType SIZE = 31;
  // lookup table for XYZ
  static constexpr uint8_t FROM_ASCII[256] =
  {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'A'     'C'             'G'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 'T'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'a'     'c'             'g'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 't'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  };

  // reverse lookup table for XYZ
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T'  // = 3
  };

  // reverse lookup table for DNA5
  static constexpr uint8_t TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0  // = 3
  };
};
struct Bits6
{
  static constexpr AlphabetSizeType SIZE = 63;
  // lookup table for XYZ
  static constexpr uint8_t FROM_ASCII[256] =
  {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'A'     'C'             'G'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 'T'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'a'     'c'             'g'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 't'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  };

  // reverse lookup table for XYZ
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T'  // = 3
  };

  // reverse lookup table for DNA5
  static constexpr uint8_t TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0  // = 3
  };
};
struct Bits7
{
  static constexpr AlphabetSizeType SIZE = 127;
  // lookup table for XYZ
  static constexpr uint8_t FROM_ASCII[256] =
  {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'A'     'C'             'G'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 'T'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'a'     'c'             'g'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 't'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  };

  // reverse lookup table for XYZ
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T'  // = 3
  };

  // reverse lookup table for DNA5
  static constexpr uint8_t TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0  // = 3
  };
};
struct Bits8
{
  static constexpr AlphabetSizeType SIZE = 255;
  // lookup table for XYZ
  static constexpr uint8_t FROM_ASCII[256] =
  {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'A'     'C'             'G'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 'T'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'a'     'c'             'g'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 't'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  };

  // reverse lookup table for XYZ
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T'  // = 3
  };

  // reverse lookup table for DNA5
  static constexpr uint8_t TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0  // = 3
  };
};



// templated test function
template<typename kmer_word_type, unsigned int kmer_size=31, class Alphabet=DNA>
void test_kmer_with_word_type_packed(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned step=1) {

  typedef MyKmer<kmer_size, Alphabet, kmer_word_type> kmer_type;

  // create (fill) Kmer
  kmer_type kmer;

  // the expected value is a 64bit kmer.  we only need the prefix corresponding to the kmer_size.
  constexpr size_t expected_shift = ((64 / bliss::AlphabetTraits<Alphabet>::getBitsPerChar() - kmer_size) * bliss::AlphabetTraits<Alphabet>::getBitsPerChar());

  unsigned char* kmer_pointer = kmer_data;
  // fill first kmer
  unsigned int offset = 0;
  offset = kmer.fillFromPackedStream(kmer_pointer, offset);

  uint64_t expected = *kmer_ex >> expected_shift;
  kmer_type kmer_ex_0(reinterpret_cast<kmer_word_type*>(&expected));
  EXPECT_EQ(kmer, kmer_ex_0) << "Kmer from stream should be equal to kmer from non-stream";


  // generate more kmers
  for (unsigned int i = step; i < nkmers; i += step)
  {
    kmer.nextFromPackedStream(kmer_pointer, offset);
    expected = kmer_ex[i] >> expected_shift;
    kmer_type kmer_ex_i(reinterpret_cast<kmer_word_type*>(&expected));
    EXPECT_EQ(kmer_ex_i, kmer) << "Kmer compare unequal for sizeof(input)="<< sizeof(unsigned char) << ", sizeof(kmer_word)=" << sizeof(kmer_word_type) << ", size=" << kmer_size << ", bits=" << bliss::AlphabetTraits<Alphabet>::getBitsPerChar() << " i = " << i;
  }
}


template<unsigned int kmer_size=31, class Alphabet=DNA>
void test_kmers_with_packed_input(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned int step=1)
{
  // test with various kmer base types
  test_kmer_with_word_type_packed<uint8_t,  kmer_size, DNA>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_packed<uint16_t, kmer_size, DNA>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_packed<uint32_t, kmer_size, DNA>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type_packed<uint64_t, kmer_size, DNA>(kmer_data, kmer_ex, nkmers, step);
}

void test_kmers_packed(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // test for bits per character: 2, 4, and 8 (no padding only!)
  test_kmers_with_packed_input<31, DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<28, DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<13, DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<4,  DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<1,  DNA>(kmer_data, kmer_ex, nkmers);

//  test_kmers_with_packed_input<10, Bits4>(kmer_data, kmer_ex, nkmers, 2);
//  test_kmers_with_packed_input<13, Bits4>(kmer_data, kmer_ex, nkmers, 2);
//
//  test_kmers_with_packed_input<7, Bits8>(kmer_data, kmer_ex, nkmers, 4);
//  test_kmers_with_packed_input<5, Bits8>(kmer_data, kmer_ex, nkmers, 4);
//
}

void test_kmers_3_packed(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // maximum in 64 bits is 21
  test_kmers_with_packed_input<21, DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<20, DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<13, DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<9,  DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<1,  DNA5>(kmer_data, kmer_ex, nkmers);
}

void test_kmers_5_packed(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // maximum in 64 bits is 12
  test_kmers_with_packed_input<12, Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<11, Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<10, Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<9,  Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<5,  Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<3,  Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_packed_input<1,  Bits5>(kmer_data, kmer_ex, nkmers);
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
   * print(",\n".join(hex((((val >> (126 - (i+3) * 2)) & 0x3) << 6) | ((val >> ((126 - (i+2) * 2)) & 0x3) << 4) | (((val >> (126 - (i+1) * 2)) & 0x3) << 2) | ((val >> (126 - i * 2)) & 0x3))[:-1] for i in range(0, 64, 4)))
   */
  uint8_t kmer_data[16] = {
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

  // test with this data

  // NO padding, therefore we can test for different input types as well

  // 8 bit input
  test_kmers_packed(kmer_data, kmer_ex, 33);

//  // 16 bit input
//  test_kmers<uint16_t>(reinterpret_cast<uint16_t*>(kmer_data), kmer_ex, 25);
//
//  // 32 bit input
//  test_kmers<uint32_t>(reinterpret_cast<uint32_t*>(kmer_data), kmer_ex, 25);
//
//  // 64 bit input
//  test_kmers<uint64_t>(reinterpret_cast<uint64_t*>(kmer_data), kmer_ex, 25);
}

/**
 * Test k-mer generation with 3 bits and thus padded input
 */
TEST(KmerGeneration, TestKmerGenerationPacked3)
{
  // test sequence: 0xabba56781234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabba56781234deadbeef01c0ffee
   * print(",\n".join([" "*24 + "0x" + hex(val >> 3*i)[-16:] for i in range(0,17)]))
   */
  uint64_t kmer_ex[17] = {0xdeadbeef01c0ffee,
                          0x9bd5b7dde0381ffd,
                          0xd37ab6fbbc0703ff,
                          0x1a6f56df7780e07f,
                          0x234deadbeef01c0f,
                          0x2469bd5b7dde0381,
                          0x048d37ab6fbbc070,
                          0xc091a6f56df7780e,
                          0x781234deadbeef01,
                          0xcf02469bd5b7dde0,
                          0x59e048d37ab6fbbc,
                          0x2b3c091a6f56df77,
                          0xa56781234deadbee,
                          0x74acf02469bd5b7d,
                          0xee959e048d37ab6f,
                          0x5dd2b3c091a6f56d,
                          0xabba56781234dead};

  // unpadded stream (bits_per_char is 3 => 1 bit padding)
  /* python:
   * ", ".join(hex((val >> i*15) & 0x7fff) for i in range(0, 9))
   */
  uint16_t kmer_data[9] = {0x7fee, 0x381, 0x7bbc, 0x756d, 0x234d, 0x4f02, 0x6e95, 0x55, 0x0};

  // test with this data
  test_kmers_3<uint16_t>(kmer_data, kmer_ex, 17);

  /* python:
   *  ", ".join(hex((val >> i*6) & 0x3f) for i in range(0, 20))
   */
  uint8_t kmer_data_8[20] = {0x2e, 0x3f, 0xf, 0x30, 0x1, 0x3c, 0x2e, 0x2f, 0x2d, 0x3a, 0xd, 0xd, 0x12, 0x20, 0x27, 0x15, 0x3a, 0x2e, 0xa, 0x0};
  // test with 8 bit (padded by 2 bits)
  test_kmers_3<uint8_t>(kmer_data_8, kmer_ex, 17);

  /* python:
   * ", ".join(hex((val >> i*30) & 0x3fffffff) for i in range(0, 5))
   */
  // 32 bits:
  uint32_t kmer_data_32[5] = {0x1c0ffee, 0x3ab6fbbc, 0x2781234d, 0x2aee95, 0x0};
  test_kmers_3<uint32_t>(kmer_data_32, kmer_ex, 17);

  /* python:
   * ", ".join(hex((val >> i*63) & 0x7fffffffffffffff) for i in range(0, 3))
   */
  // 64 bits:
  uint64_t kmer_data_64[3] = {0x5eadbeef01c0ffee, 0x15774acf02469, 0x0};
  test_kmers_3<uint64_t>(kmer_data_64, kmer_ex, 17);
}

/**
 * Test k-mer generation with 5 bits and thus padded input
 */
TEST(KmerGeneration, TestKmerGenerationPacked5)
{
  // test sequence: 0xabba56781234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabba56781234deadbeef01c0ffee
   * print(",\n".join([" "*24 + "0x" + hex(val >> 5*i)[-16:] for i in range(0,11)]))
   */
  uint64_t kmer_ex[11] = {0xdeadbeef01c0ffee,
                          0xa6f56df7780e07ff,
                          0x8d37ab6fbbc0703f,
                          0x2469bd5b7dde0381,
                          0x81234deadbeef01c,
                          0x3c091a6f56df7780,
                          0x59e048d37ab6fbbc,
                          0x4acf02469bd5b7dd,
                          0xba56781234deadbe,
                          0x5dd2b3c091a6f56d,
                          0x2aee959e048d37ab};

  // unpadded stream (bits_per_char is 5 => 1 bit padding)
  /* python:
   * ", ".join(hex((val >> i*15) & 0x7fff) for i in range(0, 9))
   */
  uint16_t kmer_data[9] = {0x7fee, 0x381, 0x7bbc, 0x756d, 0x234d, 0x4f02, 0x6e95, 0x55, 0x0};

  // test with this data
  test_kmers_5<uint16_t>(kmer_data, kmer_ex, 11);

  /* python:
   *  ", ".join(hex((val >> i*5) & 0x1f) for i in range(0, 24))
   */
  uint8_t kmer_data_8[24] = {0xe, 0x1f, 0x1f, 0x1, 0x1c, 0x0, 0x1c, 0x1d, 0x1e, 0xd, 0xb, 0x1d, 0xd, 0x1a, 0x8, 0x2, 0x18, 0x13, 0x15, 0x14, 0x1b, 0x15, 0x2, 0x0};
  // test with 8 bit (padded by 3 bits)
  test_kmers_5<uint8_t>(kmer_data_8, kmer_ex, 11);

  /* python:
   * ", ".join(hex((val >> i*30) & 0x3fffffff) for i in range(0, 5))
   */
  // 32 bits (2 bits padding)
  uint32_t kmer_data_32[5] = {0x1c0ffee, 0x3ab6fbbc, 0x2781234d, 0x2aee95, 0x0};
  test_kmers_5<uint32_t>(kmer_data_32, kmer_ex, 11);

  /* python:
   * ", ".join(hex((val >> i*60) & 0x0fffffffffffffff) for i in range(0, 3))
   */
  // 64 bits (4 bits padding)
  uint64_t kmer_data_64[3] = {0xeadbeef01c0ffee, 0xabba56781234d, 0x0};
  test_kmers_5<uint64_t>(kmer_data_64, kmer_ex, 11);
}



// templated test function
template<typename kmer_word_type, unsigned int kmer_size=31, class Alphabet=DNA>
void test_kmer_with_word_type(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned step=1) {

  typedef MyKmer<kmer_size, Alphabet, kmer_word_type> kmer_type;

  // create (fill) Kmer
  kmer_type kmer;

  // the expected value is a 64bit kmer.  we only need the prefix corresponding to the kmer_size.
  constexpr size_t expected_shift = ((64 / bliss::AlphabetTraits<Alphabet>::getBitsPerChar() - kmer_size) * bliss::AlphabetTraits<Alphabet>::getBitsPerChar());
  //printf("expected_shift: %lu\n", expected_shift);

  unsigned char* kmer_pointer = kmer_data;
  //printf("kmer pointer: %d\n", *kmer_pointer);
  // fill first kmer
  //unsigned int offset = 0;
  kmer.fillFromChars(kmer_pointer);
  uint64_t expected = *kmer_ex >> expected_shift;
  kmer_type kmer_ex_0(reinterpret_cast<kmer_word_type*>(&expected));
  //printf("iter i = %d, expected: %016lX, actual %016lX\n", 0, kmer_ex_0.getPrefix64(), kmer.getPrefix64());

  EXPECT_EQ(kmer, kmer_ex_0) << "Kmer from stream should be equal to kmer from non-stream";


  // generate more kmers
  for (unsigned int i = step; i < nkmers; i += step)
  {
    kmer.nextFromChar(*kmer_pointer); ++kmer_pointer;
    expected = kmer_ex[i] >> expected_shift;
    kmer_type kmer_ex_i(reinterpret_cast<kmer_word_type*>(&expected));

    //printf("iter i = %d, expected: %016lX, actual %016lX\n", i, kmer_ex_i.getPrefix64(), kmer.getPrefix64());

    EXPECT_EQ(kmer_ex_i, kmer) << "Kmer compare unequal for sizeof(input)="<< sizeof(unsigned char) << ", sizeof(kmer_word)=" << sizeof(kmer_word_type) << ", size=" << kmer_size << ", bits=" << bliss::AlphabetTraits<Alphabet>::getBitsPerChar() << " i = " << i;
  }
}


template<unsigned int kmer_size=31, class Alphabet=DNA>
void test_kmers_with_char_input(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers, unsigned int step=1)
{
  // test with various kmer base types
  test_kmer_with_word_type<uint8_t,  kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type<uint16_t, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type<uint32_t, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
  test_kmer_with_word_type<uint64_t, kmer_size, Alphabet>(kmer_data, kmer_ex, nkmers, step);
}

void test_kmers(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // test for bits per character: 2, 4, and 8 (no padding only!)
  test_kmers_with_char_input<31, DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<28, DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<13, DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<4,  DNA>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<1,  DNA>(kmer_data, kmer_ex, nkmers);
}

void test_kmers_3(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // maximum in 64 bits is 21
  test_kmers_with_char_input<21, DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<20, DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<13, DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<9,  DNA5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<1,  DNA5>(kmer_data, kmer_ex, nkmers);
}

void test_kmers_5(unsigned char* kmer_data, uint64_t* kmer_ex, unsigned int nkmers)
{
  // maximum in 64 bits is 12
  test_kmers_with_char_input<12, Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<11, Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<10, Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<9,  Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<5,  Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<3,  Bits5>(kmer_data, kmer_ex, nkmers);
  test_kmers_with_char_input<1,  Bits5>(kmer_data, kmer_ex, nkmers);
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
  test_kmers(kmer_data, kmer_ex, 33);

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
  test_kmers_3(kmer_data_8, kmer_ex, 22);

}

/**
 * Test k-mer generation with 5 bits and thus padded input
 */
TEST(KmerGeneration, TestKmerGenerationChar5)
{
  // test sequence: 0xabbacafebabe1234deadbeef01c0ffee

  // expected kmers:
  // generated by the python commands (thank you python for integrated bigints)
  /*
   * val = 0xabbacafebabe1234deadbeef01c0ffee
   * print(",\n".join([" "*24 + hex(val >> (128 - 60 - 5*i) & 0x0fffffffffffffff)[:-1] for i in range(0,14)]))
   */
  uint64_t kmer_ex[14] = {
                          0xabbacafebabe123,
                          0x77595fd757c2469,
                          0xeb2bfaeaf848d37,
                          0x657f5d5f091a6f5,
                          0xafebabe1234dead,
                          0xfd757c2469bd5b7,
                          0xaeaf848d37ab6fb,
                          0xd5f091a6f56df77,
                          0xbe1234deadbeef0,
                          0xc2469bd5b7dde03,
                          0x48d37ab6fbbc070,
                          0x1a6f56df7780e07,
                          0x4deadbeef01c0ff,
                          0xbd5b7dde0381ffd
};


  /* python:
   * print(",\n".join([" "*24 + hex(val >> (128 - 5 - 5*i) & 0x1f)[:-1] for i in range(0,25)]))
   */
  uint8_t kmer_data_8[25] = {
                             0x15,
                             0xe,
                             0x1d,
                             0xc,
                             0x15,
                             0x1f,
                             0x15,
                             0x1a,
                             0x17,
                             0x18,
                             0x9,
                             0x3,
                             0x9,
                             0x17,
                             0x15,
                             0xd,
                             0x17,
                             0x1b,
                             0x17,
                             0x10,
                             0x3,
                             0x10,
                             0x7,
                             0x1f,
                             0x1d
  };
  // test with 8 bit (padded by 3 bits)
  test_kmers_5(kmer_data_8, kmer_ex, 14);

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

  MyKmer<41, DNA, uint16_t> kmer(kmer_val);
  MyKmer<41, DNA, uint16_t> kmer_s(kmer_val_s4);
  MyKmer<41, DNA, uint16_t> kmer_g(kmer_val_g3);

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
  uint16_t kmer_val[] = {0xffee, 0x1c0, 0xbeef, 0xdead, 0x1234, 0x5678, 0xabba};
  uint16_t kmer_ex[] = {0xaeea, 0x2d95, 0x1c84, 0x7ab7, 0xfbbe, 0x340, 0xbbff};
  uint16_t kmer_ex_3[] = {0x2faa, 0x2795, 0x34a4, 0xdabd, 0x3ebe, 0x2311, 0x6bff};
  uint16_t kmer_ex_4[] = {0xabba, 0x8765, 0x4321, 0xdaed, 0xfeeb, 0xc10, 0xeeff};
  uint16_t kmer_ex_5[] = {0xd375, 0xb13a, 0xba40, 0xd5f5, 0xe77c, 0x8780, 0x1dff};
  uint16_t kmer_ex_7[] = {0xb755, 0xcf2, 0xa644, 0xd6bd, 0x1777, 0x18ee, 0xddfc};


  /* test for bits_per_char = 2 */
  MyKmer<7*8, DNA, uint16_t> kmer_in(kmer_val);
  MyKmer<7*8, DNA, uint16_t> kmer_ex_rev(kmer_ex);
  // get the reverse
  MyKmer<7*8, DNA, uint16_t> kmer_rev = kmer_in.reversed_kmer();
  EXPECT_EQ(kmer_ex_rev, kmer_rev);

  /* test for bits_per_char = 3 */
  MyKmer<37, DNA5, uint16_t> kmer3_in(kmer_val);
  MyKmer<37, DNA5, uint16_t> kmer3_ex_rev(kmer_ex_3);
  // get the reverse
  MyKmer<37, DNA5, uint16_t> kmer3_rev = kmer3_in.reversed_kmer();
  EXPECT_EQ(kmer3_ex_rev, kmer3_rev);

  /* test for bits_per_char = 4 */
  MyKmer<28, Bits4, uint16_t> kmer4_in(kmer_val);
  MyKmer<28, Bits4, uint16_t> kmer4_ex_rev(kmer_ex_4);
  // get the reverse
  MyKmer<28, Bits4, uint16_t> kmer4_rev = kmer4_in.reversed_kmer();
  EXPECT_EQ(kmer4_ex_rev, kmer4_rev);

  /* test for bits_per_char = 5 */
  MyKmer<22, Bits5, uint16_t> kmer5_in(kmer_val);
  MyKmer<22, Bits5, uint16_t> kmer5_ex_rev(kmer_ex_5);
  // get the reverse
  MyKmer<22, Bits5, uint16_t> kmer5_rev = kmer5_in.reversed_kmer();
  EXPECT_EQ(kmer5_ex_rev, kmer5_rev);

  /* test for bits_per_char = 7 */
  MyKmer<16, Bits7, uint16_t> kmer7_in(kmer_val);
  MyKmer<16, Bits7, uint16_t> kmer7_ex_rev(kmer_ex_7);
  // get the reverse
  MyKmer<16, Bits7, uint16_t> kmer7_rev = kmer7_in.reversed_kmer();
  EXPECT_EQ(kmer7_ex_rev, kmer7_rev);
  // TODO: tests for smaller sizes and bigger bits_per_char

}
