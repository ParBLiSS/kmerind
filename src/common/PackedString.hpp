/**
 * @file    PackedString.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares the PackedStringImpl and PackedString classes.
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_PACKEDSTRING_H
#define BLISS_COMMON_PACKEDSTRING_H


// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <vector>
#include <iterator> // for std::begin() and std::end()

// own includes
#include <common/base_types.hpp>
#include <common/bit_ops.hpp>


template <BitSizeType bitsPerChar, typename CharType = CharType, typename WordType = WordType>
class PackedStringImpl
{
public:
  typedef size_t size_type;
  typedef size_type index_type;
  typedef CharType value_type;
  // base container
  // TODO consider using raw memory
  std::vector<WordType> packedString;

  size_type nChars;
  static constexpr BitSizeType charsPerWord = (sizeof(WordType) * 8) / bitsPerChar;
  size_type nWords;
  static constexpr BitSizeType BITS_PER_CHAR = bitsPerChar;

  inline void packIntoWord(WordType& word, const CharType c, const BitSizeType nBits, const BitSizeType pos)
  {
    word |= static_cast<WordType>(getCharBitMask(nBits) & c) << (nBits * pos);
  }

  PackedStringImpl()
  {

  }

  template <typename container_type>
  PackedStringImpl(const container_type& container)
  {
    this->packChars(container);
  }

  template <typename IteratorType>
  PackedStringImpl(const IteratorType begin, const IteratorType end)
  {
    this->packChars(begin, end);
  }

  template <typename container_type>
  void packChars(const container_type& container)
  {
    this->packChars(std::begin(container), std::end(container));
  }

  // TODO could work with dynamic arrays without a randomaccessIterator
  template <typename IteratorType>
  void packChars(const IteratorType begin, const IteratorType end)
  {
    // TODO: this requires a RandomAccessIterator, packing could though also work from a simple forward iterator using
    // the dynamic array property of std::vector
    this->nChars = end - begin;

    // get the needed number of words (emulating ceil(nChars/charsPerWord))
    nWords = nChars / charsPerWord + (nChars % charsPerWord != 0 ? 1 : 0);

    // initialize all values to 0
    this->packedString = std::vector<WordType>(nWords, 0);

    // get iterators
    auto it = begin;
    auto out_it = this->packedString.begin();

    for (unsigned int i = 0; i < nChars;)
    {
      // pack all chars into the current word
      packIntoWord(*out_it, *it, bitsPerChar, i % charsPerWord);

      // update iterators
      ++it;
      ++i;
      if (i % charsPerWord == 0) ++out_it;
    }
  }

  inline CharType unpackChar(const WordType& word, const BitSizeType char_idx) const
  {
    WordType mask = getWordBitMask(bitsPerChar, char_idx * bitsPerChar);
    return static_cast<CharType>((word & mask) >> (char_idx * bitsPerChar));
  }

  // TODO: is std::vector the best container in this case?
  template <typename OutputIterator>
  void unpackSequence(OutputIterator outIter) const
  {
    size_type j = 0;
    for (WordType word : packedString)
    {
      for (unsigned int i = 0; i < charsPerWord; ++i)
      {
        if (j++ == this->size()) return;
        *(outIter++) = unpackChar(word, i);
      }
    }
  }

public:

  const value_type operator[](index_type idx) const
  {
    size_type word_idx = idx / charsPerWord;
    BitSizeType char_idx = idx % charsPerWord;
    return unpackChar(packedString[word_idx], char_idx);
  }

  size_type size() const
  {
    return this->nChars;
  }
};

// TODO consider using alphabet traits, instead of bits_per_char_needed
template <typename T>
using PackedString = PackedStringImpl<bits_per_char_needed(T::SIZE)>;

#endif // BLISS_COMMON_PACKEDSTRING_H
