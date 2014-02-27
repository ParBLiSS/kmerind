/**
 * @file    Kmer.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Implements the Kmer class.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_KMER_H
#define BLISS_COMMON_KMER_H

// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <iterator>
#include <type_traits>

// own includes
#include <common/base_types.hpp>
#include <common/bit_ops.hpp>
#include <common/padding.hpp>


namespace bliss
{


template <unsigned int KMER_SIZE, unsigned int BITS_PER_CHAR, typename word_type=WordType>
class Kmer
{
private:
  /// The size of the Kmer, i.e. the number of characters
  static constexpr unsigned int size = KMER_SIZE;
  /// The number of bits of each character
  static constexpr unsigned int bitsPerChar = BITS_PER_CHAR;
  /// The total number of bits
  static constexpr unsigned int nBits = size * bitsPerChar;

  /// The padding traits of an unpadded stream (i.e. the padding for the very
  /// last word)
  typedef UnpaddedStreamTraits<word_type, nBits> bitstream;

  /// The number of stored words inside the Kmer
  static constexpr unsigned int nWords = bitstream::nWords;

  /* last character offsets (and whether or not it is split accord storage
   * words)
   */
  /// Offset to the very last character in the k-mer
  static constexpr unsigned int lastCharOffset = nBits - bitsPerChar;
  /// Whether or not the last character is split across storage words
  static constexpr bool lastCharIsSplit = (lastCharOffset < (nWords-1)*sizeof(word_type)*8);
  /// The offset to the very last character by word boundary
  static constexpr unsigned int lastCharWordOffset = lastCharOffset % (sizeof(word_type)*8);
  /// In case the last character is split: the number of bits of the last
  /// character in the previous from last word.
  static constexpr unsigned int leftSplitSize = sizeof(word_type)*8 - lastCharWordOffset;

  /// The actual storage of the k-mer
  word_type data[nWords];

public:

  Kmer()
  {
    // make this a valid Kmer, other bits are left uninitialized
    do_sanitize();
  }

  template<typename InputIterator>
  Kmer(InputIterator begin)
  {
    static_assert(std::is_same<
        typename std::iterator_traits<InputIterator>::value_type,
        word_type>::value,
        "Input iterator must have same value type as the Kmer storage");

    // copy all the data into this kmer
    word_type* out = data;
    for (unsigned int i = 0; i < nWords; ++i)
    {
      *(out++) = *(begin++);
    }

    // set unused bits to zero to make this a valid kmer
    do_sanitize();
  }

  /*
   * TODO:
   *  - fill and nextKmer each for unpadded streams
   *  - reverse + complement functions
   *    (complement is function of the alphabet)
   *  - XOR operator (!?)
   *  - hash function (!?)
   */

  template <typename InputIterator>
  // TODO: add option for bit offset in input sequence?
  unsigned int fillFromPaddedStream(InputIterator& begin)
  {
    // remove padding copy to own data structure
    typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
    const unsigned int paddingBits = PaddingTraits<input_word_type, bitsPerChar>::padding_bits;
    removePadding(begin, data, size*bitsPerChar, paddingBits);

    // set unused bits to 0
    do_sanitize();

    // the bit offset in the input sequence
    // TODO: do this inside the removePadding function!?
    unsigned int offset = (size*bitsPerChar) % PaddingTraits<input_word_type, bitsPerChar>::data_bits;
    std::advance(begin, (size*bitsPerChar) / PaddingTraits<input_word_type, bitsPerChar>::data_bits);

    // return the offset
    return offset;
  }

  template <typename InputIterator>
  void nextFromPaddedStream(InputIterator& begin, unsigned int& offset)
  {
    typedef typename std::iterator_traits<InputIterator>::value_type input_type;
    // shift the kmer by the size of one character
    do_right_shift(bitsPerChar);

    // iterate to next word if the current one is done
    while (offset >= PaddingTraits<input_type, bitsPerChar>::data_bits)
    {
      ++begin;
      offset -= PaddingTraits<input_type, bitsPerChar>::data_bits;
    }

    // get the bigger type of the input_type and the kmer word_type
    typedef typename std::conditional<sizeof(input_type) < sizeof(word_type),
                                      word_type, input_type>::type bigger_type;

    // get the next input as the bigger word type
    bigger_type curWord = static_cast<bigger_type>(*begin);
    // shift offset between input offset and the offset for the last character
    // in the kmer
    int shift_by = offset - lastCharWordOffset;
    if (shift_by >= 0)
    {
      // positive shift: right shift
      // (shift_by == 0 is handled in this case as well, because this is the
      //  faster case, and for = 0 it doesn't make a difference which one is
      //  used)
      data[nWords - 1 - lastCharIsSplit] = static_cast<word_type>(curWord >> shift_by);
    }
    else
    {
      // negative shift: left shift
      data[nWords - 1 - lastCharIsSplit] |= static_cast<word_type>(curWord << -shift_by);
    }

    // in case the last character is split across words in the kmer
    // representation: set last bits to last word
    // NOTE: this branch is removed by the compiler (result known at compile time)
    if (lastCharIsSplit)
    {
      data[nWords - 1] = static_cast<word_type>(curWord >> (offset + leftSplitSize));
    }

    // set unused bits to 0
    do_sanitize();

    // increase offset
    offset += bitsPerChar;
  }

  bool operator==(const Kmer& rhs) const
  {
    return std::equal(data, data+nWords, rhs.data);
  }

  bool operator!=(const Kmer& rhs) const
  {
    return !(this->operator==(rhs));
  }

  // for debug purposes
  // TODO: support different alphabets
  std::string toString()
  {
    // 1.) unpacking
    // 2.) back-translating
  }

protected:

  inline void do_sanitize()
  {
    // TODO use templated helper struct for <0> template specialization
    data[nWords-1] &= getBitMask<word_type>(bitstream::invPadBits);
  }


  // TODO template specialization for nWords = 1 (just use base type shift)
  // TODO implement more efficient version doing fixed left shift by BITS_PER_CHAR
  inline void do_left_shift(size_t shift)
  {
    // inspired by STL bitset implementation
    // TODO: replace by std::div
    const size_t word_shift = shift / (sizeof(word_type)*8);
    const size_t offset = shift % (sizeof(word_type)*8);

    if (offset == 0)
    {
      // no bit shifting, just shift words around
      for (size_t i = nWords - 1; i >= word_shift; --i)
      {
        data[i] = data[i - word_shift];
      }
    }
    else
    {
      const size_t inv_offset = sizeof(word_type)*8 - offset;
      for (size_t i = nWords - 1; i > word_shift; --i)
      {
        data[i] = ((data[i - word_shift] << offset) | (data[i - word_shift - 1] >> inv_offset));
      }
      data[word_shift] = data[0] << offset;
    }
    // set all others to 0
    std::fill(data, data+word_shift, static_cast<word_type>(0));
  }

  // TODO template specialization for nWords = 1 (just use base type shift)
  // TODO implement more efficient version doing fixed left shift by BITS_PER_CHAR
  inline void do_right_shift(size_t shift)
  {
    // inspired by STL bitset implementation
    // TODO replace by std::div
    const size_t word_shift = shift / (sizeof(word_type)*8);
    const size_t offset = shift % (sizeof(word_type)*8);

    if (offset == 0)
    {
      // no bit shifting, just shift words around
      for (size_t i = 0; i < nWords - word_shift; ++i)
      {
        data[i] = data[i + word_shift];
      }
    }
    else
    {
      const size_t inv_offset = sizeof(word_type)*8 - offset;
      for (size_t i = 0; i < nWords - word_shift - 1; ++i)
      {
        data[i] = ((data[i + word_shift] >> offset) | (data[i + word_shift + 1] << inv_offset));
      }
      data[nWords - word_shift - 1] = data[nWords-1] >> offset;
    }
    // set all others to 0
    std::fill(data + (nWords - word_shift), data + nWords, static_cast<word_type>(0));
  }
};

} // namespace bliss

#endif // BLISS_COMMON_KMER_H