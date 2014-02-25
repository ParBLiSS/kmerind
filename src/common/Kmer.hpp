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
  static constexpr unsigned int size = KMER_SIZE;
  static constexpr unsigned int bitsPerChar = BITS_PER_CHAR;
  static constexpr unsigned int nBits = size * bitsPerChar;
  static constexpr unsigned int nBytes = intCeil<unsigned int>(nBits, 8);
  static constexpr unsigned int nWords = intCeil<unsigned int>(nBytes, sizeof(word_type));
  // padding internal to Kmer (only in last storage word)
  static constexpr unsigned int padBits = nWords*sizeof(word_type)*8 - nBits;
  static constexpr unsigned int invPadBits = sizeof(word_type)*8 - padBits;
  static constexpr unsigned int padBytes = nWords*sizeof(word_type) - nBytes;

  // last character offsets (and whether or not it is split accord storage words)
  static constexpr unsigned int lastCharOffset = nBits - bitsPerChar;
  static constexpr bool lastCharIsSplit = (lastCharOffset < (nWords-1)*sizeof(word_type)*8);
  static constexpr unsigned int lastCharWordOffset = lastCharOffset % (sizeof(word_type)*8);
  static constexpr unsigned int leftSplitSize = sizeof(word_type)*8 - lastCharWordOffset;


  // storage:
  word_type data[nWords];

public:

  Kmer()
  {
    // make this a valid Kmer
    do_sanitize();
  }

  template<typename InputIterator>
  Kmer(InputIterator begin)
  {
    static_assert(std::is_same<typename std::iterator_traits<InputIterator>::value_type, word_type>::value, "Input iterator must have same value type as the Kmer storage");
    word_type* out = data;
    for (unsigned int i = 0; i < nWords; ++i)
    {
      *(out++) = *(begin++);
    }
    do_sanitize();
  }

  inline void do_sanitize()
  {
    // TODO use templated helper struct for <0> template specialization
    data[nWords-1] &= getBitMask<word_type>(invPadBits);
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


  template <typename InputIterator>
  Kmer& fillFromPaddedStream(InputIterator begin)
  {
    // remove padding copy to own data structure
    // TODO somehow save current offset for nextKmerFromPaddedStream() ?
    typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
    const unsigned int paddingBits = PaddingTraits<input_word_type, bitsPerChar>::padding_bits;
    removePadding(begin, data, size*bitsPerChar, paddingBits);
    // set unused bits to 0
    do_sanitize();

    return (*this);
  }


  template <typename InputIterator>
  Kmer& nextKmerFromPaddedStream(InputIterator& begin, unsigned int& offset)
  {
    typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
    // 1.) right shift this kmer by BITS_PER_CHAR
    do_right_shift(bitsPerChar);

    // iterate to next word if the current one is done
    while (offset >= PaddingTraits<input_word_type, bitsPerChar>::data_bits)
    {
      ++begin;
      offset -= PaddingTraits<input_word_type, bitsPerChar>::data_bits;
    }

    typedef typename std::conditional<sizeof(input_word_type) < sizeof(word_type), word_type, input_word_type>::type bigger_word_type;

    bigger_word_type curWord = static_cast<bigger_word_type>(*begin);
    int shift_by = offset - lastCharWordOffset;
    if (shift_by > 0)
    {
      data[nWords - 1 - lastCharIsSplit] = static_cast<word_type>(curWord >> shift_by);
    }
    else
    {
      data[nWords - 1 - lastCharIsSplit] |= static_cast<word_type>(curWord << -shift_by);
    }
    // NOTE: this branch is removed by the compiler (result known at compile time)
    if (lastCharIsSplit)
    {
      data[nWords - 1] = static_cast<word_type>(curWord >> (offset + leftSplitSize));
    }
    // set unused bits to 0
    do_sanitize();

    // increase offset
    offset += bitsPerChar;

    return (*this);
  }

  bool operator==(const Kmer& rhs) const
  {
    return std::equal(data, data+nWords, rhs.data);
  }

  bool operator!=(const Kmer& rhs) const
  {
    return !(this->operator==(rhs));
  }

};

} // namespace bliss

#endif // BLISS_COMMON_KMER_H