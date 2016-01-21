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
 * @file    packed_string.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares the PackedStringImpl and PackedString classes.
 */
#ifndef BLISS_COMMON_PACKEDSTRING_H
#define BLISS_COMMON_PACKEDSTRING_H

// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <vector>
#include <iterator> // for std::begin() and std::end()

// own includes
#include "common/base_types.hpp"
#include "common/bit_ops.hpp"

namespace bliss
{
  namespace common
  {
  /**
   * @brief A class implementing a PackedString.
   *
   * The class is templated by at least the bit size into which characters are
   * to be packed. E.g. for an alphabet of size 4 (e.g. DNA) this parameter
   * would be set to 2 bits.
   *
   * The parameters are template parameters:
   *
   * @tparam BITS_PER_CHAR  The number of bits needed for each packed character.
   * @tparam CharType       The underlying type of the unpacked characters.
   *                        Default = char.
   * @tparam WordType       The underlying type used for storage. This is usually
   *                        set to the biggest machine word (i.e. uint64_t for
   *                        a 64 bit system).
   */
  template <BitSizeType BITS_PER_CHAR, typename CharType = CharType, typename WordType = WordType>
  class PackedStringImpl
  {
  public:
    /// An unsigned integral type representing the size of the PackedString
    typedef size_t size_type;
  
    /// An unsigned integral type representing an index into the PackeString
    typedef size_type index_type;
  
    /// The type of the unpacked characters, returned by the [] operator and
    /// the unpackSequence() function.
    typedef CharType value_type;
  
    /**
     * @brief   Constructs a new PackedString of size
     *          `container.end() - container.begin()` and packs the given
     *          sequence into the new PackedString.
     *
     * Any linear/contiguous STL (or custom) container can be used. This includes
     * std::string, std::basic_string<DNA>, std::vector, std::array, etc.
     *
     * @param container[in]   A STL container that implements the `begin()` and
     *                        `end()` methods.
     * @tparam ContainerType  A STL container type.
     */
    template <typename ContainerType>
    PackedStringImpl(const ContainerType& container)
    {
      this->packChars(container);
    }
  
    /**
     * @brief   Constructs a new PackedString of size `end - begin` and packs the
     *          given sequence into this new PackedString.
     *
     * @param begin   An InputIterator to the first position of the sequence to
     *                be packed and stored in this PackedString sequence.
     * @param end     An InputIterator to the final position (one after the last
     *                element) of the sequence to be packed and stored in this
     *                PackedString.
     *
     * @tparam IteratorType   InputIterator type for the `begin()` and `end()`
     *                        iterators.
     */
    template <typename IteratorType>
    PackedStringImpl(const IteratorType begin, const IteratorType end)
    {
      this->packChars(begin, end);
    }
  
    /**
     * @brief Unpacks the PackedString into the given OutputIterator.
     *
     * The OutputIterator will be increased `size()` times. Please make sure
     * the according container is large enough to hold `size()` elements.
     *
     * @param outIter   An OutputIterator to the initial position of the
     *                  destination sequence.
     * @returns         An iterator to the end of the destination sequence.
     *                  This iterator will point towards one element past
     *                  the last written element.
     * @tparam  OutputIterator  A OutputIterator type that can be used in
     *                          sequential output operations.
     */
    template <typename OutputIterator>
    OutputIterator unpackSequence(OutputIterator outIter) const
    {
      // counter of elements written out
      size_type j = 0;
      // loop through all words and all chars per word
      for (WordType word : packedString)
      {
        for (unsigned int i = 0; i < charsPerWord; ++i)
        {
          // quit iterating in case all elements have been unpacked
          if (j++ == this->size())
            return outIter;
          // unpack the next element and increase output iterator
          *(outIter++) = unpackChar(word, i);
        }
      }
      return outIter;
    }
  
    /**
     * @brief Returns the value of the packed char with position `idx`.
     *
     * @param idx   The index of the packed character to be returned.
     * @return      The packed character with position `idx`.
     */
    const value_type operator[](index_type idx) const
    {
      // get the index of the word
      size_type word_idx = idx / charsPerWord;
      // get the index of the char within the word
      BitSizeType char_idx = idx % charsPerWord;
      // return the unpacked char
      return unpackChar(packedString[word_idx], char_idx);
    }
  
    /**
     * @brief   Returns the size (i.e. the number of characters stored in the
     *          PackedString) of the PackedString.
     *
     * @return  The size of the PackedString.
     */
    size_type size() const
    {
      return this->nChars;
    }
  
  private:
    /// The base container for the packed sequence.
    std::vector<WordType> packedString;
  
    /// The number of characters in the packed sequence. This is the total size
    /// of the sequence.
    size_type nChars;
  
    /// The number of bits that one char is packed into
    static constexpr BitSizeType bitsPerChar = BITS_PER_CHAR;
  
    /// The number of characters that are packed into one storage word.
    static constexpr BitSizeType charsPerWord = (sizeof(WordType) * 8) / bitsPerChar;
  
    /// The number of storage words. This is equivalent to the size of the base
    /// container.
    size_type nWords;
  
    // disallow the default constructor (by making it private)
    PackedStringImpl() {}
  
  
    /**
     * @brief Packs the given character into the given word at the given position.
     *
     * This function packs a character into the given storage word. For this
     * the character is cut off to BITS_PER_CHAR bits, where only the least
     * significant bits are used. The offset given by the parameter `pos` starts
     * at the least significant bit of the storage word. E.g. for the value of
     * pos = 0, the character will stored into the `BITS_PER_CHAR` least
     * significant bits of the storage word. For each other pos = i, the
     * character will first be shifted by `i*BITS_PER_CHAR` prior to storage
     * inside the storage word.
     *
     * Only the bits [i*BITS_PER_CHAR, (i+1)*BITS_PER_CHAR - 1] are changed inside
     * of `word`.
     *
     * @param word[in|out]  The given character is packed into this word.
     * @param c[in]         The character to be packed.
     * @param pos[in]       The position inside the word of where the character
     *                      will be packed to.
     */
    inline void packIntoWord(WordType& word, const CharType c, const BitSizeType pos)
    {
      word |= static_cast<WordType>(getCharBitMask(bitsPerChar) & c) << (bitsPerChar * pos);
    }
  
    /**
     * @brief Packs a sequence given by the STL container into this PackedString.
     *
     * @param container   A STL container type that implements the iterator
     *                    functions `begin()` and `end()`.
     * @tparam ContainerType  A STL container type.
     */
    template <typename container_type>
    void packChars(const container_type& container)
    {
      this->packChars(std::begin(container), std::end(container));
    }
  
    /**
     * @brief Packs a sequence given by the `begin` and `end` iterators into
     *        this PackedString.
     * @param begin   An InputIterator to the first position of the sequence to
     *                be packed and stored in this PackedString sequence.
     * @param end     An InputIterator to the final position (one after the last
     *                element) of the sequence to be packed and stored in this
     *                PackedString.
     * @tparam IteratorType   InputIterator type for the `begin()` and `end()`
     *                        iterators.
     */
    template <typename IteratorType>
    void packChars(const IteratorType begin, const IteratorType end)
    {
      // TODO: this requires a RandomAccessIterator
      // packing could though also work from a simple forward iterator using
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
        packIntoWord(*out_it, *it, i % charsPerWord);
  
        // update iterators
        ++it;
        ++i;
        if (i % charsPerWord == 0) ++out_it;
      }
    }
  
    /**
     * @brief Unpacks a single character from the given word at the given index.
     *
     * @param word[in]      The storage word from which the character is to be
     *                      unpacked.
     * @param char_idx[in]  The index to the character to be unpacked and
     *                      returned.
     */
    inline CharType unpackChar(const WordType& word, const BitSizeType char_idx) const
    {
      WordType mask = getWordBitMask(bitsPerChar, char_idx * bitsPerChar);
      return static_cast<CharType>((word & mask) >> (char_idx * bitsPerChar));
    }
  };
  
  /**
   * @brief A packed string representation for the given alphabet T.
   *
   * This is a wrapper type around the PackedStringImpl. This wrapper
   * is templated by the alphabet (e.g. DNA, DNA5, AA, etc.) while the
   * actual implementation PackedStringImpl is templated by the number
   * of bits used for storage.
   *
   * Parameters are template parameters:
   *
   * @tparam T    The alphabet type. At the least, this class must implement
   *              the static ::SIZE attribute, which returns the size of the
   *              alphabet.
   */
  // TODO consider using alphabet traits, instead of bits_per_char_needed
  template <typename T>
  using PackedString = PackedStringImpl<ceilLog2(T::SIZE)>;
  
  } // namespace common
} // namespace bliss

#endif // BLISS_COMMON_PACKEDSTRING_H
