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
 * @file    kmer.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the Kmer data type.
 *
 */
#ifndef BLISS_COMMON_KMER_H
#define BLISS_COMMON_KMER_H

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

#define KMER_INLINE inline

namespace bliss
{

  namespace common
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
  template <unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE=WordType>
  class Kmer
  {

      static_assert(::std::is_integral<WORD_TYPE>::value && !::std::is_signed<WORD_TYPE>::value, "Kmer only accepts primitive unsigned integral type for internal storage.");

      friend std::string bliss::utils::KmerUtils::toASCIIString< Kmer >(const Kmer & kmer);
  
      template<unsigned int K, typename ALPHA, typename WT>
      friend std::ostream& operator<<(std::ostream& ost, const Kmer<K, ALPHA, WT> & kmer);



   public:
      /// The size of the Kmer, i.e. the number of characters
    static constexpr unsigned int size = KMER_SIZE;
    /// The number of bits of each character
    static constexpr unsigned int bitsPerChar = bliss::common::AlphabetTraits<ALPHABET>::getBitsPerChar();
  
    typedef WORD_TYPE KmerWordType;
    typedef ALPHABET KmerAlphabet;

    /// The total number of bits
    static constexpr unsigned int nBits = size * bitsPerChar;

   private:
  
    /// The padding traits of an unpadded stream (i.e. the padding for the very
    /// last word)
    typedef UnpaddedStreamTraits<WORD_TYPE, nBits> bitstream;
  
   public:
    /// The number of stored words inside the Kmer
    static constexpr unsigned int nWords = bitstream::nWords;
    static constexpr unsigned int nBytes = bitstream::nBytes;
  
   private:
  
    /*
     * last character offsets (and whether or not it is split accord storage
     * words)
     */
    /// Offset to the very last character in the k-mer
    static constexpr unsigned int lastCharOffset = nBits - bitsPerChar;
    /// Whether or not the last character is split across storage words
    static constexpr bool lastCharIsSplit = (lastCharOffset < (nWords-1)*sizeof(WORD_TYPE)*8);
    /// The offset to the very last character by word boundary
    static constexpr unsigned int lastCharWordOffset = lastCharOffset % (sizeof(WORD_TYPE)*8);
    /// In case the last character is split: the number of bits of the last
    /// character in the previous from last word.
    static constexpr unsigned int leftSplitSize = sizeof(WORD_TYPE)*8 - lastCharWordOffset;

    using MACH_WORD_TYPE = size_t;

    /// The actual storage of the k-mer
    WORD_TYPE data[nWords]; // alignas is not compatible with mxx datatype stuff - causes size of std::pair to increase greatly because alignment of elements is based on most resitrictive element.
  
   public:

    /**
     * @brief   The default constructor, creates an uninitialized k-mer.
     *
     * The k-mer's data is not initalized, but still represents a valid k-mer,
     * as the unused bits are set to 0.
     */
    explicit Kmer(bool clear = true)
    {
      // make this a valid Kmer, other bits are left uninitialized
      // do_sanitize();
  
      // set all bits to zero (no more 'may be used uninitialized' warning)
      if (clear) do_clear();
    }
  
    /// copy constructor  - want result of iterator to be copiable here so no "explicit"
    Kmer(Kmer const& other) {
      memcpy(this->data, other.data, nWords * sizeof(WORD_TYPE));
    };
  
    /// copy assignment operator
    Kmer& operator=(Kmer const &other) {

      memcpy(this->data, other.data, nWords * sizeof(WORD_TYPE));
      return *this;
    }


    /*
     * TODO:
     *  - fill and nextKmer each for unpadded streams
     *  - reverse + complement functions
     *    (complement is function of the alphabet)
     *  - hash function (!?)
     */

    /**
     * @brief get the internal data for direct access.
     *        useful for hash function and others.
     * @return pointer to internal data
     */
    WORD_TYPE const * getData() const  {
      return data;
    }

    WORD_TYPE* getData() {
    	return data;
    }

    WORD_TYPE (&getDataRef())[nWords] {
      return data;
    }

    /**
     * FIXME: update documentation for `size => (size-1)`
     * @brief   Fills this k-mer from a packed and padded input sequence.
     * @note    one example of padded input sequence is where the character values are in
     *          a contiguous range (e.g. DNA 0..3), stored in a char, so 2 bits are used, and 6 bits are padding.
     *
     *          another example of packed and padded input sequence is where the character values are in
     *          a contiguous range DNA5, 0..4), stored in 2 in a char, and 3 bits are used, 2 bits are padding.
     *
     *          packed word has ordering such that higher order bits have characters later in sequence.
     *            e.g. sequence ACGT may have a packed word representation with TGCA
     *          i.e. prefix at the LSB side.
     *
     * The k-mer's data is filled from the given input sequence using
     * `KMER_SIZE * BITS_PER_CHAR` bits.
     *
     * @tparam InputIterator  An input iterator type that is at least a forward
     *                        iterator.
     * @param begin[in|out]   An interator reference pointing to the beginning of
     *                        the packed and padded sequence used to fill the
     *                        k-mer's data structure. When this function returns,
     *                        this iterator will point to the next element of the
     *                        sequence to be read. Together with the bit offset
     *                        returned by this function, this defines where to
     *                        continue to generate k-mers from the input sequence
     *                        via the `nextFromPaddedStream()` function.
     * @param stop_on_last    Whether to stop the iterator on the last read
     *                        element, rather than stopping on the next element
     *                        after the one that was last read. This influences
     *                        both the given iterator `begin` and the returned
     *                        offset.
     * @returns               The bit offset of the current iterator position
     *                        to the bits that have to be read in the next
     *                        iteration.
     */
    template <typename InputIterator, typename offset_t>
    // TODO: add option for bit offset in input sequence?
    unsigned int fillFromPackedStream(InputIterator& begin, offset_t& offset, bool stop_on_last = false)
    {
  
      typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
  //    const unsigned int paddingBits = PackingTraits<input_word_type, bitsPerChar>::padding_bits;
  //    removePadding(begin, data, size*bitsPerChar, paddingBits);
  
      //static_assert(sizeof(input_word_type) == 1, "only support byte as PackStream data type");
  
      const unsigned int input_data_bits = PackingTraits<input_word_type, bitsPerChar>::data_bits;
  
      // iterate to next word if the current one is done
      while (offset >= input_data_bits)
      {
        ++begin;
        offset -= input_data_bits;
      }
  
      // clear k-mer
      this->do_clear();
  
      int bitPos = bitstream::nBits - bitsPerChar;
  
      // add to lsb iteratively, one packed word at a time
      for (size_t i = 0; i < KMER_SIZE; ++i, bitPos -= bitsPerChar)
      {
        setBitsAtPos(*begin >> offset, bitPos, bitsPerChar);
  
        // don't move iterator during last iteration if that option is set
        if (i == (KMER_SIZE - 1) && stop_on_last) continue;
  
        // increase offset
        offset += bitsPerChar;
        if (offset >= input_data_bits)
        {
          ++begin;
          offset -= input_data_bits;
        }
      }
  
      // set unused bits to 0
      do_sanitize();
  
  //    // the bit offset in the input sequence
  //    // TODO: do this inside the removePadding function!?
  //    unsigned int total_bits = stop_on_last ? (size-1)*bitsPerChar : size*bitsPerChar;
  //    unsigned int offset = total_bits % PackingTraits<input_word_type, bitsPerChar>::data_bits;
  //    std::advance(begin, total_bits / PackingTraits<input_word_type, bitsPerChar>::data_bits);
  
      // return the offset
      return offset;
    }
  
  
    /**
     * FIXME: update documentation for `size => (size-1)`
     * @brief   Fills this k-mer from a packed and padded input sequence.
     * @note    one example of padded input sequence is where the character values are in
     *          a contiguous range (e.g. DNA 0..3), stored in a char, so 2 bits are used, and 6 bits are padding.
     *
     *          another example of packed and padded input sequence is where the character values are in
     *          a contiguous range DNA5, 0..4), stored in 2 in a char, and 3 bits are used, 2 bits are padding.
     *
     *          packed word has ordering such that higher order bits have characters later in sequence.
     *            e.g. sequence ACGT may have a packed word representation with TGCA
     *          i.e. prefix at the LSB side.
     *
     * The k-mer's data is filled from the given input sequence using
     * `KMER_SIZE * BITS_PER_CHAR` bits.
     *
     * @tparam InputIterator  An input iterator type that is at least a forward
     *                        iterator.
     * @param begin[in|out]   An interator reference pointing to the beginning of
     *                        the packed and padded sequence used to fill the
     *                        k-mer's data structure. When this function returns,
     *                        this iterator will point to the next element of the
     *                        sequence to be read. Together with the bit offset
     *                        returned by this function, this defines where to
     *                        continue to generate k-mers from the input sequence
     *                        via the `nextFromPaddedStream()` function.
     * @param stop_on_last    Whether to stop the iterator on the last read
     *                        element, rather than stopping on the next element
     *                        after the one that was last read. This influences
     *                        both the given iterator `begin` and the returned
     *                        offset.
     * @returns               The bit offset of the current iterator position
     *                        to the bits that have to be read in the next
     *                        iteration.
     */
    template <typename InputIterator, typename offset_t>
    // TODO: add option for bit offset in input sequence?
    unsigned int fillReverseFromPackedStream(InputIterator& begin, offset_t& offset, bool stop_on_last = false)
    {
  
      typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
  //    const unsigned int paddingBits = PackingTraits<input_word_type, bitsPerChar>::padding_bits;
  //    removePadding(begin, data, size*bitsPerChar, paddingBits);
  
      //static_assert(sizeof(input_word_type) == 1, "only support byte as PackStream data type");
  
      const unsigned int input_data_bits = PackingTraits<input_word_type, bitsPerChar>::data_bits;
  
      // iterate to next word if the current one is done
      while (offset >= input_data_bits)
      {
        ++begin;
        offset -= input_data_bits;
      }
  
      // clear k-mer
      this->do_clear();
  
      // TODO:  can copy directly from stream into kmer, IF data types are same.
      int bitPos = 0;
      // add to lsb iteratively, one packed word at a time
      for (unsigned int i = 0; i < KMER_SIZE; ++i, bitPos += bitsPerChar)
      {
        setBitsAtPos(*begin >> offset, bitPos, bitsPerChar);
  
        // don't move iterator during last iteration if that option is set
        if (i == (KMER_SIZE - 1) && stop_on_last) continue;
  
        // increase offset
        offset += bitsPerChar;
        if (offset >= input_data_bits)
        {
          ++begin;
          offset -= input_data_bits;
        }
      }
  
      // set unused bits to 0
      do_sanitize();
  
  //    // the bit offset in the input sequence
  //    // TODO: do this inside the removePadding function!?
  //    unsigned int total_bits = stop_on_last ? (size-1)*bitsPerChar : size*bitsPerChar;
  //    unsigned int offset = total_bits % PackingTraits<input_word_type, bitsPerChar>::data_bits;
  //    std::advance(begin, total_bits / PackingTraits<input_word_type, bitsPerChar>::data_bits);
  
      // return the offset
      return offset;
    }
  
  /**
   * Note that iterator has value domain consistent with the valid values in the alphabet
   * @param begin
   * @param stop_on_last
   */
    template <typename InputIterator>
    void fillFromChars(InputIterator& begin, bool stop_on_last = false)
    {
  
      // clear k-mer
      this->do_clear();
  
      // directly put the char where it needs to go.
  
      int bitPos = bitstream::nBits - bitsPerChar;
      for (size_t i = 0; i < KMER_SIZE; ++i, bitPos -= bitsPerChar) {
  
        setBitsAtPos(*begin, bitPos, bitsPerChar);
  
        // don't move iterator during last iteration if that option is set
        if (i == (KMER_SIZE - 1) && stop_on_last) continue;
  
        ++begin;
      }
  
    }
  
    /**
     * @brief fill kmer from InputIterator, where each element contains 1 single character.
     * @note  input iterator has values that are valid for the alphabet.
     * @param begin
     * @param stop_on_last
     */
    template <typename InputIterator>
    void fillReverseFromChars(InputIterator& begin, bool stop_on_last = false)
    {
  
      // clear k-mer
      this->do_clear();
  
      // directly put the char where it needs to go.
  
      int bitPos = 0;
      for (unsigned int i = 0; i < KMER_SIZE; ++i, bitPos += bitsPerChar) {
  
        setBitsAtPos(*begin, bitPos, bitsPerChar);
  
        // don't move iterator during last iteration if that option is set
        if (i == (KMER_SIZE - 1) && stop_on_last) continue;
  
        ++begin;
      }
  
    }
  
    /**
     * @brief   Generates the next k-mer from the given sequence using a sliding
     *          window approach.
     *
     * Given the packed and padded input sequence via the `begin` iterator
     * and the current bit offset, this function will read the next
     * `BITS_PER_CHAR` bits, left shift the current k-mer value by that number
     * of bits and puts the newly read bits into the `BITS_PER_CHAR` least
     * significant bits of the k-mer value. Both parameters are passed by
     * reference and internally updated to the next position to read from.
     * Therefore, the user of this function just needs to keep calling this
     * function without updating the parameters to generate all k-mers for a
     * given packed and padded input sequence.
     *
     * Note that each element in the InputIterator contains an integral number of
     * packed characters starting from the lsb position.
     *
     * @tparam InputIterator  An input iterator type that is at least a forward
     *                        iterator.
     * @param begin           An iterator pointing to the next position of the
     *                        input sequence to be read.
     * @param offset          The bit offset of where to start reading inside
     *                        the currently pointed to value of the input
     *                        sequence.
     */
    template <typename InputIterator, typename offset_t>
    void nextFromPackedStream(InputIterator& begin, offset_t& offset)
    {
      typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
  //  //     shift the kmer by the size of one character
  //  //     TODO: replace this by a call that does exactly bitsPerChar right shift
  //  //     (better compiler optimization)
  //    do_left_shift(bitsPerChar);
  
      //static_assert(sizeof(input_type) == 1, "only support byte as PackStream data type");
  
      const unsigned int input_data_bits = PackingTraits<input_word_type, bitsPerChar>::data_bits;
  
      // iterate to next word if the current one is done
      while (offset >= input_data_bits)
      {
        ++begin;
        offset -= input_data_bits;
      }
  
      nextFromWordInternal(*begin >> offset, bitsPerChar);
  
      // set unused bits to 0
      do_sanitize();
  
      // increase offset
      offset += bitsPerChar;
      if (offset >= input_data_bits)
      {
        ++begin;
        offset -= input_data_bits;
      }
    }
  
    /**
     * @brief   Generates the next k-mer from the given sequence using a sliding
     *          window approach.
     *
     * Given the packed and padded input sequence via the `begin` iterator
     * and the current bit offset, this function will read the next
     * `BITS_PER_CHAR` bits, left shift the current k-mer value by that number
     * of bits and puts the newly read bits into the `BITS_PER_CHAR` least
     * significant bits of the k-mer value. Both parameters are passed by
     * reference and internally updated to the next position to read from.
     * Therefore, the user of this function just needs to keep calling this
     * function without updating the parameters to generate all k-mers for a
     * given packed and padded input sequence.
     *
     * Note that each element in the InputIterator contains an integral number of
     * packed characters starting from the lsb position.
     *
     * @tparam InputIterator  An input iterator type that is at least a forward
     *                        iterator.
     * @param begin           An iterator pointing to the next position of the
     *                        input sequence to be read.
     * @param offset          The bit offset of where to start reading inside
     *                        the currently pointed to value of the input
     *                        sequence.
     */
    template <typename InputIterator, typename offset_t>
    void nextReverseFromPackedStream(InputIterator& begin, offset_t& offset)
    {
      typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
  //  //     shift the kmer by the size of one character
  //  //     TODO: replace this by a call that does exactly bitsPerChar right shift
  //  //     (better compiler optimization)
  //    do_left_shift(bitsPerChar);
  
      //static_assert(sizeof(input_type) == 1, "only support byte as PackStream data type");
  
      const unsigned int input_data_bits = PackingTraits<input_word_type, bitsPerChar>::data_bits;
  
      // iterate to next word if the current one is done
      while (offset >= input_data_bits)
      {
        ++begin;
        offset -= input_data_bits;
      }
  
      nextReverseFromWordInternal(*begin >> offset, bitsPerChar);
  
      // set unused bits to 0
      do_sanitize();
  
      // increase offset
      offset += bitsPerChar;
      if (offset >= input_data_bits)
      {
        ++begin;
        offset -= input_data_bits;
      }
    }
  
  
    // TODO:  need a version of nextFromChar that takes an iterator and increments from the iterator.
  
  
    /**
     * @brief Creates the next k-mer by shifting and adding the given character.
     *
     * Creates the next k-mer by the sliding window. This first left shifts the
     * current k-mer and then adds the given character to the least significant
     * bits of this k-mer.
     *
     * @note  input iterator has values that are valid for the alphabet.
     * Note: this is changing this k-mer in-place.
     *
     * @param c The character to be added to the k-mer.
     */
    void nextFromChar(unsigned char c)
    {
      // shift the kmer by the size of one character
      // TODO: replace this by a call that does exactly bitsPerChar right shift
      // (better compiler optimization)
  
      nextFromWordInternal(c, bitsPerChar);

      // clean up
      do_sanitize();
    }
  
    // TODO: generating the reverse .
  
    /**
     * @brief Creates the next reverse k-mer by shifting and adding the given character.
     *
     * Creates the next k-mer by the sliding window. This first right shifts the
     * current k-mer and then adds the given character to the most significant
     * bits of this k-mer.
     *
     * @note  input iterator has values that are valid for the alphabet.
     *
     * Note: this is changing this k-mer in-place.
     *
     * @param c The character to be added to the k-mer.
     */
    void nextReverseFromChar(unsigned char c)
    {
      // shift the kmer by the size of one character
      // TODO: replace this by a call that does exactly bitsPerChar right shift
      // (better compiler optimization)
  
      nextReverseFromWordInternal(c, bitsPerChar);
  
      // clean up  - not needed - internally shifting to the right.
      do_sanitize();
    }
  
  
    /* equality comparison operators */
  
    /**
     * @brief Compares this k-mer with the given k-mer for equality.
     *
     * Returns whether the data of the k-mers is identical.
     *
     * @returns   `true` if the data of both k-mers is identical, `false`
     *            otherwise.
     */
    KMER_INLINE bool operator==(const Kmer& rhs) const
    {
      return (memcmp(data, rhs.data, nBytes) == 0);
    }
  
    /**
     * @brief Compares this k-mer with the given k-mer for in-equality.
     *
     * Returns whether the data of the k-mers is different.
     *
     * @returns   `false` if the data of both k-mers is identical, `true`
     *            otherwise.
     */
    KMER_INLINE bool operator!=(const Kmer& rhs) const
    {
      return !(this->operator==(rhs));
    }
  
    /* ordered comparison operators */
  
    /**
     * @brief Returns whether this k-mer compares smaller than the given k-mer.
     * @note  kmers have "least recently seen" character at most significant bit.
     *        and lexicographic comparison is order from least recently seen (prefix to suffix)
     *        so comparison needs to go from MSB to LSB.
     *
     * @returns `True` if this k-mer compares smaller than the given k-mer,
     *          `False` otherwise.
     */
    KMER_INLINE bool operator<(const Kmer& rhs) const
    {
      auto first = this->data + nWords - 1;
      auto second = rhs.data + nWords - 1;

      for (size_t i = 0; i < nWords; ++i, --first, --second) {
        if (*first != *second) return (*first < *second);  // if equal, keep comparing. else decide.
      }
      return false;  // all equal
//      std::pair<const WORD_TYPE*, const WORD_TYPE*> unequal = std::mismatch(this->data, this->data + nWords, rhs.data);
//      if (unequal.first == this->data + nWords)
//      {
//        // all elements are equal
//        return false;
//      }
//      else
//      {
//        // the comparison of the first unequal element will determine the result
//        return *(unequal.first) < *(unequal.second);
//      }
//      return memcmp(data, rhs.data, nWords * sizeof(WORD_TYPE)) < 0;
    }
  
    /**
     * @brief Returns whether this k-mer compares smaller or equal than the given
     *        k-mer.
     *
     * @returns `True` if this k-mer compares smaller than or equal to the given
     *          k-mer, `False` otherwise.
     */
    KMER_INLINE bool operator<=(const Kmer& rhs) const
    {
      return !(rhs < *this);
    }
  
    /* symmetric comparison operators */
  
    /**
     * @brief Returns whether this k-mer compares greater or equal than the given
     *        k-mer.
     *
     * @returns `True` if this k-mer compares greater than or equal to the given
     *          k-mer, `False` otherwise.
     */
    KMER_INLINE bool operator>=(const Kmer& rhs) const
    {
      return !(*this < rhs);
    }
  
    /**
     * @brief Returns whether this k-mer compares greater than the given k-mer.
     *
     * @returns `True` if this k-mer compares greater than the given k-mer,
     *          `False` otherwise.
     */
    KMER_INLINE bool operator>(const Kmer& rhs) const
    {
      return rhs < *this;
    }
  
    /* bit operators */
  
    /**
     * @brief XOR
     */
    KMER_INLINE Kmer& operator^=(const Kmer& rhs)
    {
      for (uint32_t i = 0; i < nWords; ++i) this->data[i] ^= rhs.data[i];
//      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_xor<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief XOR
     */
    KMER_INLINE Kmer operator^(const Kmer& rhs) const
    {
      Kmer result = *this;
      result ^= rhs;
      return result;
    }
  
    /**
     * @brief AND
     */
    KMER_INLINE Kmer& operator&=(const Kmer& rhs)
    {
      for (uint32_t i = 0; i < nWords; ++i) this->data[i] &= rhs.data[i];
//      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_and<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief AND
     */
    KMER_INLINE Kmer operator&(const Kmer& rhs) const
    {
      Kmer result = *this;
      result &= rhs;
      return result;
    }
  
    /**
     * @brief OR
     */
    KMER_INLINE Kmer& operator|=(const Kmer& rhs)
    {
      for (uint32_t i = 0; i < nWords; ++i) this->data[i] |= rhs.data[i];
//      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_or<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief OR
     */
    KMER_INLINE Kmer operator|(const Kmer& rhs) const
    {
      Kmer result = *this;
      result |= rhs;
      return result;
    }
  
  
    /**
     * @brief Shifts the k-mer left by the given number of CHARACTERS.
     *
     * @note  This shifts by the number of characters (which is a larger shift
     *        then bitwise).
     */
    KMER_INLINE Kmer& operator<<=(const std::size_t shift_by)
    {
      // shift the given number of _characters_!
      this->do_left_shift(shift_by * bitsPerChar);
      this->do_sanitize();
      return *this;
    }
  
    /**
     * @brief Shifts the k-mer left by the given number of CHARACTERS.
     *
     * @note  This shifts by the number of characters (which is a larger shift
     *        then bitwise).
     */
    KMER_INLINE Kmer operator<<(const std::size_t shift_by)
    {
      Kmer result = *this;
      result <<= shift_by;
      return result;
    }
  
    /**
     * @brief Shifts the kmer left by the given number of bits.
     * @param shift_by
     * @return
     */
    KMER_INLINE Kmer left_shift_bits(const std::size_t shift_by) {
      this->do_left_shift(shift_by);
      this->do_sanitize();
      return *this;
    }

    /**
     * @brief Shifts the k-mer right by the given number of CHARACTERS.
     *
     * @note  This shifts by the number of characters (which is a larger shift
     *        then bitwise).
     */
    KMER_INLINE Kmer& operator>>=(const std::size_t shift_by)
    {
      this->do_sanitize();  // right shift so need to sanitize upper bits first.

      // shift the given number of **characters** !
      this->do_right_shift(shift_by * bitsPerChar);
      return *this;
    }
  
  
    /**
     * @brief Shifts the k-mer right by the given number of CHARACTERS.
     *
     * @note  This shifts by the number of characters (which is a larger shift
     *        then bitwise).
     */
    KMER_INLINE Kmer operator>>(const std::size_t shift_by)
    {
      Kmer result = *this;
      result >>= shift_by;
      return result;
    }
  
    /**
     * @brief Shifts the kmer right by the given number of bits.
     * @param shift_by
     * @return
     */
    KMER_INLINE Kmer right_shift_bits(const std::size_t shift_by) {
      this->do_right_shift(shift_by);
      this->do_sanitize();
      return *this;
    }

  
    /**
     * @brief Returns a reversed k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    KMER_INLINE Kmer reverse() const
    {
      // only needs to know the bitsPerChar info.
      return do_reverse(*this);
    }

    /**
     * @brief Returns a reverse complement of a k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    KMER_INLINE Kmer reverse_complement() const
    {
      return do_reverse_complement<ALPHABET>(*this);
    }
  

    // for debug purposes
    /**
     * @brief Returns a string representation of this k-mer.
     *
     * Usage: for debugging and testing.
     *
     * @returns   A std::string representing this k-mer.
     */
    std::string toString() const
    {
      /* return the hex representation of the data array values */
      std::stringstream ss;
      ss << "k-mer of size " << size << " alphabet ";
      ss << typeid(ALPHABET).name();
      ss << ": MSB [";
      for (int64_t i = nWords - 1; i >=0; --i)
      {
        ss << "0x" << std::hex << data[i] << " ";
      }
      ss << "]";
      return ss.str();
    }
  
    // for debug purposes
    /**
     * @brief Returns a string representation of this k-mer.  Order is 5' to 3' (correspond to MSB to LSB)
     *
     * Usage: for debugging and testing.
     *
     * @returns   A std::string representing this k-mer.
     */
    std::string toAlphabetString() const
    {
      /* return the char representation of the data array values */
      std::stringstream result;
      Kmer cpy(*this);

      std::stack<size_t> elementsInReverse;
      size_t forBitMask = (1 << bitsPerChar) - 1;


      for (unsigned int i = 0; i < size; ++i)
      {
        elementsInReverse.push(static_cast<size_t>(forBitMask & cpy.data[0]));
        cpy.do_right_shift(bitsPerChar);
      }

      //Pop stack elements to string stream
      while(!elementsInReverse.empty())
      {
        result << elementsInReverse.top() << " ";
        elementsInReverse.pop();
      }

      return result.str();
    }



    /// get 64 bit prefix.  for hashing.
    uint64_t getPrefix(const unsigned int NumBits = 64) const {
  
      if (bitstream::padBits + NumBits <= bitstream::bitsPerWord) {
        // word contains the NumBits bits and the padding bits.
        // need to shift to the right to get rid of extra bits
        return static_cast<uint64_t>(data[nWords - 1] >>
                                     (bitstream::bitsPerWord - NumBits - bitstream::padBits));
  
      } else {
        // padBits + NumBits is larger than a word.
          // more than 1 word, so shift and get suffix
          Kmer temp = *this;
          if (bitstream::nBits > NumBits) {
            temp.do_right_shift(bitstream::nBits - NumBits);
          }
          return temp.getSuffix(NumBits);
      }
  
    }
  
  
  
    uint64_t getInfix(const unsigned int NumBits = 64, const unsigned int offsetFromMSB = 0) const {
      if (bitstream::nBits - NumBits > offsetFromMSB ) {
  
        // more than 64 bits to the right of offsetFromMSB.  so need to shift to right.
        Kmer temp = *this;
        temp.do_right_shift(bitstream::nBits - offsetFromMSB - NumBits);
        return temp.getSuffix(NumBits);
      } else if (bitstream::nBits - NumBits == offsetFromMSB) {
        // exactly 64 bits left.  just use it.
        return getSuffix(NumBits);
      } else {
  
        // less than 64 bits left at offset from MSB.  don't need to shift to right.
        // but do need to clear leading portion.
        return getSuffix(NumBits) & getLeastSignificantBitsMask<uint64_t>(bitstream::nBits - offsetFromMSB);
  
      }
  
    }
  
    uint64_t getSuffix(const unsigned int NumBits = 64) const {
      if (sizeof(WORD_TYPE) * 8 >= NumBits || bitstream::nWords == 1) {
        // kmer composes of 1+ words that are larger or equal to 64 bits, or a single word.
        return static_cast<uint64_t>(data[0]) & getLeastSignificantBitsMask<uint64_t>(NumBits);
  
  
      } else {
        // kmer has multiple small words. compose it.
        const size_t nwords = static_cast<size_t>(bitstream::nWords) <= ((NumBits + sizeof(WORD_TYPE) * 8 - 1) / (sizeof(WORD_TYPE) * 8) ) ?
            static_cast<size_t>(bitstream::nWords) :
            ((NumBits + sizeof(WORD_TYPE) * 8 - 1) / (sizeof(WORD_TYPE) * 8) );
        uint64_t result = 0;
        for (int64_t i = nwords - 1; i >= 0; --i) {
  
          // this is to avoid GCC compiler warning.  even though the case
          // sizeof (WORD_TYPE) >= sizeof(uint64_t) is already caught earlier, compiler
          // will still check this line and cause a warning to be thrown.
          // so we artificially insert a conditional that caps WORD_TYPE size to 7
          // of word type will never go above 7 (actually 4) at runtime (caught by branch earlier)
          // and most of this code will be optimized out as well during compilation.
          result <<= ((sizeof(WORD_TYPE) > 7 ? 7 : sizeof(WORD_TYPE)) * 8);
          result |= data[i];
        }
        return result;
      }
    }
  
  
  
  protected:
  
    template<typename WType>
    KMER_INLINE void setBitsAtPos(WType w, int bitPos, unsigned int numBits) {
      // get the value masked.
      WORD_TYPE charVal = static_cast<WORD_TYPE>(w) & getLeastSignificantBitsMask<WORD_TYPE>(numBits);
  
      // determine which word in kmer it needs to go to
      unsigned int wordId = bitPos / bitstream::bitsPerWord;
      unsigned int offsetInWord = bitPos % bitstream::bitsPerWord;  // offset is where the LSB of the char will sit, in bit coordinate.
  
      if (offsetInWord >= (bitstream::bitsPerWord - numBits)) {
        // if split between words, deal with it.
        data[wordId] |= charVal << offsetInWord;   // the lower bits of the charVal
  
        // the number of lowerbits consumed is (bitstream::bitsPerWord - offsetInWord)
        // so right shift those many places and what remains goes into the next word.
        if (wordId < nWords - 1) data[wordId + 1] |= charVal >> (bitstream::bitsPerWord - offsetInWord);
  
  
      } else {
        // else insert into the specific word.
        data[wordId] |= charVal << offsetInWord;
      }
  
    }
  
    /**
     * @brief internal method to add one more character to the kmer at the LSB side
     * @param c     character to add.
     */
    template <typename WType>
    KMER_INLINE void nextFromWordInternal(WType w, unsigned int shift)
    {
      // left shift k-mer
      // TODO: replace by single shift operation
      do_left_shift(shift);
  
      // add character to least significant end (requires least shifting)
      *data |= static_cast<WORD_TYPE>(w) &
          getLeastSignificantBitsMask<WORD_TYPE>(shift);
    }

    /**
     * @brief internal method to add one more character to the kmer at the MSB side
     * @param c     character to add.
     */
    template <typename WType>
    KMER_INLINE void nextReverseFromWordInternal(WType w, unsigned int shift)
    {
      // left shift k-mer
      // TODO: replace by single shift operation
      do_right_shift(shift);
  
      // add character to least significant end (requires least shifting)
      data[nWords - 1] |= (static_cast<WORD_TYPE>(w) &
          getLeastSignificantBitsMask<WORD_TYPE>(shift)) << (bitstream::invPadBits - shift);
    }
  
    /**
     * @brief Sets all unused bits of the underlying k-mer data to 0.
     * @details  highest order bits in highest number element are 0.
     */
    KMER_INLINE void do_sanitize()
    {
      // TODO use templated helper struct for <0> template specialization
      data[nWords-1] &= getLeastSignificantBitsMask<WORD_TYPE>(bitstream::invPadBits);
    }

    /**
     * @brief Sets all bits to zero.
     */
    KMER_INLINE void do_clear()
    {
      //std::fill(data, data + nWords, static_cast<WORD_TYPE>(0));
      memset(data, 0, nWords * sizeof(WORD_TYPE));
    }
  
    /**
     * @brief Performs a left shift by `shift` bits on the k-mer data.
     *
     * @param shift   The number of bits to shift by.
     */
    // TODO template specialization for nWords = 1 (just use base type shift)
    // TODO implement more efficient version doing fixed left shift by BITS_PER_CHAR
    KMER_INLINE void do_left_shift(size_t const & shift)
    {
      // inspired by STL bitset implementation
      const int64_t word_shift = shift / (sizeof(WORD_TYPE)*8);
      const size_t offset = shift % (sizeof(WORD_TYPE)*8);
  
      // all shifted away.
      if (word_shift >= nWords) {
        do_clear();
        return;
      }
  
      if (offset == 0)
      {
        // no bit shifting, just shift words around.  do in reverse so don't overwrite.
        for (int64_t i = nWords - 1; i >= word_shift; --i)
        {
          data[i] = data[i - word_shift];
        }
      }
      else
      {
        const size_t inv_offset = sizeof(WORD_TYPE)*8 - offset;
        WORD_TYPE t = data[nWords - 1 - word_shift];
        WORD_TYPE t1;

        for (int64_t i = nWords - 1, j = nWords - word_shift - 2; i > word_shift; --i, --j)
        {
          t1 = data[j];
          data[i] = ((t << offset) | (t1 >> inv_offset));
          t = t1;
        }
        data[word_shift] = data[0] << offset;
      }
  
      // set all others to 0
//      std::fill(data, data+word_shift, static_cast<WORD_TYPE>(0));
      memset(data, 0, word_shift * sizeof(WORD_TYPE));
    }
  
    /**
     * @brief Performs a right shift by `shift` bits on the k-mer data.
     *
     * @param shift   The number of bits to shift by.
     */
    // TODO template specialization for nWords = 1 (just use base type shift)
    // TODO implement more efficient version doing fixed right shift by BITS_PER_CHAR
    KMER_INLINE void do_right_shift(size_t const & shift)
    {
      // inspired by STL bitset implementation
      const size_t word_shift = shift / (sizeof(WORD_TYPE)*8);
      const size_t offset = shift % (sizeof(WORD_TYPE)*8);
  
      // all shifted away.
      if (word_shift >= nWords) {
        do_clear();
        return;
      }
  
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
        const size_t inv_offset = sizeof(WORD_TYPE)*8 - offset;
        WORD_TYPE t = data[word_shift];
        WORD_TYPE t1;
        for (size_t i = 1; i < nWords - word_shift; ++i)
        {
          t1 = data[word_shift + i];
          data[i-1] = ((t >> offset) | (t1 << inv_offset));
          t = t1;
        }
        data[nWords - word_shift - 1] = t >> offset;
      }
      // set all others to 0
//      std::fill(data + (nWords - word_shift), data + nWords, static_cast<WORD_TYPE>(0));
      memset(data + (nWords - word_shift), 0, word_shift * sizeof(WORD_TYPE));
    }
  




    /// reverse the kmer.  only parameter of ALPHABET that matters here is bitsPerChar.
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, DNA6>::value ||
                                  ::std::is_same<A, RNA6>::value, int>::type = 0>
    KMER_INLINE Kmer do_reverse(Kmer const & src) const
    {
      Kmer result(false);

//      //const uint8_t* in = reinterpret_cast<const uint8_t*>(src.data);
//      uint8_t* out = reinterpret_cast<uint8_t*>(result.data);
//
//      // NOTE: AVX2 is not faster.
//      ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_AVX2>(out, reinterpret_cast<const uint8_t*>(src.data), nBytes);
//
//      // TODO: more efficient right shift.
//      result.do_right_shift(nBytes * 8 - nBits);
      if (nWords * sizeof(WORD_TYPE) >= 16)  // choose based on performance.  AVX is faster >256 bit, below it SWAR is faster.
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_AVX2>(result.data, src.data);
      else if (nWords * sizeof(WORD_TYPE) >= 12)
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SSSE3>(result.data, src.data);
      else
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SWAR>(result.data, src.data);

      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return result;
    }

    /// reverse the kmer.  only parameter of ALPHABET that matters here is bitsPerChar.
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA6>::value ||
                                    ::std::is_same<A, RNA6>::value), int>::type = 0>
    KMER_INLINE Kmer do_reverse(Kmer const & src) const
    {
      Kmer result(false);

//      //const uint8_t* in = reinterpret_cast<const uint8_t*>(src.data);
//      uint8_t* out = reinterpret_cast<uint8_t*>(result.data);
//
//      // NOTE: AVX2 is not faster.
//      ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_AVX2>(out, reinterpret_cast<const uint8_t*>(src.data), nBytes);
//
//      // TODO: more efficient right shift.
//      result.do_right_shift(nBytes * 8 - nBits);
      if (nWords * sizeof(WORD_TYPE) > 32)  // choose based on performance.  AVX is faster >256 bit (and fallback to SSSE3 is okay), below 256 SWAR is faster.
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_AVX2>(result.data, src.data);
//      else if (nWords * sizeof(WORD_TYPE) > 16)
//        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SSSE3>(result.data, src.data);
      else
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SWAR>(result.data, src.data);

      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return result;
    }


    /// reverse complement of kmer.  specialzied for DNA/RNA, where the complement is the bitwise negation.
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value, int>::type = 0>
    KMER_INLINE Kmer do_reverse_complement(Kmer const& src) const
    {
      // repeat code from do_reverse to avoid separate do_sanitize.

      // DNA and RNA complement is via negation.
      // working with words here since we want to speed up negation.
      Kmer result(false);

      //const uint8_t* in = reinterpret_cast<const uint8_t*>(src.data);
//      uint8_t* out = reinterpret_cast<uint8_t*>(result.data);

      // NOTE: SWAR is faster for size <= 256bit, after which AVX2 (and SSSE3) is faster.
      if (nWords * sizeof(WORD_TYPE) > 32)
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_AVX2>(result.data, src.data);
//      else if (nWords * sizeof(WORD_TYPE) > 16)
//        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SSSE3>(result.data, src.data);
      else
        ::bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SWAR>(result.data, src.data);
      //for (uint32_t i = 0; i < nWords; ++i) result.data[i] = ~result.data[i];
      ::bliss::utils::bit_ops::negate(result.data, result.data, nWords);

      // TODO: more efficient right shift.
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return result;  // negate by xor with 0
    }

    /// reverse complement of kmer.  specialzied for DNA5/DNA6/RNA5/RNA6/DNA16, where the complement is the bitwise reverse.
    template <typename A = ALPHABET,
        typename ::std::enable_if< ::std::is_same<A, DNA6>::value ||
                                  ::std::is_same<A, RNA6>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    KMER_INLINE Kmer do_reverse_complement(Kmer const& src) const
    {
      // DNA6, RNA6, and DNA16 use 1 bit reverse.   DNA5 and RNA5 are aliased to DNA6 and RNA6

      Kmer result(false);
//      result.data[nWords - 1] = 0;

      //const uint8_t* in = reinterpret_cast<const uint8_t*>(src.data);
//      uint8_t* out = reinterpret_cast<uint8_t*>(result.data);
      if (nWords * sizeof(WORD_TYPE) >= 48)
      ::bliss::utils::bit_ops::reverse<1, ::bliss::utils::bit_ops::BIT_REV_AVX2>(result.data, src.data);
      else if (nWords * sizeof(WORD_TYPE) >= 32)
      ::bliss::utils::bit_ops::reverse<1, ::bliss::utils::bit_ops::BIT_REV_SSSE3>(result.data, src.data);
      else
      ::bliss::utils::bit_ops::reverse<1, ::bliss::utils::bit_ops::BIT_REV_SWAR>(result.data, src.data);


      // TODO: MORE EFFICIENT SHIFT
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);
      return result;
    }

    /// other alphabet types - reverse complement via serial implementation and applies ALPHABET's TO_COMPLEMENT lookup.
    /// if word size is multiple of bits per char (so basically, bits per char is power of 2.
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA6>::value ||
                                    ::std::is_same<A, RNA6>::value ||
                                    ::std::is_same<A, DNA16>::value) &&
                                     ((sizeof(WORD_TYPE) * 8) % bitsPerChar != 0), int>::type = 0>
    KMER_INLINE Kmer do_reverse_complement(Kmer const & src) const
    {
      static_assert(!(::std::is_same<A, DNA>::value ||
          ::std::is_same<A, RNA>::value ||
          ::std::is_same<A, DNA6>::value ||
          ::std::is_same<A, RNA6>::value ||
          ::std::is_same<A, DNA16>::value) &&
           ((sizeof(WORD_TYPE) * 8) % bitsPerChar != 0), "do reverse complement is not defined for alphabet with size != 2, 3, 4 and word type not a multiple of bits Per char.");
    }

//    template <typename A = ALPHABET,
//        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
//                                    ::std::is_same<A, RNA>::value ||
//                                    ::std::is_same<A, DNA6>::value ||
//                                    ::std::is_same<A, RNA6>::value ||
//                                    ::std::is_same<A, DNA16>::value) &&
//                                     ((sizeof(WORD_TYPE) * 8) % bitsPerChar == 0), int>::type = 0>
//    KMER_INLINE Kmer do_reverse_complement(Kmer const & src) const
//    {
//
//      /* Linear (inefficient) reverse:  because we don't have a specialization that supports TO_COMPLEMENT.  do word by word..*/
//
//      // get temporary copy of this
//      Kmer comp, result(false);
//
//      // get lower most bits from the temp copy and push them into the lower bits
//      // of this
//      WORD_TYPE tmp, tmp2;
//
//      ::std::cout << ::std::hex;
//
//      // complement all, then reverse.
//      for (unsigned int i = 0; i < nWords; ++i)
//      {
//        tmp = src.data[i];
//
//
//        for (unsigned int j = 0; j < (sizeof(WORD_TYPE) * 8); j += bitsPerChar) {
////          ::std::cout << "i " << i << " j " << j << " tmp: " << tmp << " ";
//          tmp2 = tmp & (getLeastSignificantBitsMask<WORD_TYPE>(bitsPerChar) << j);  // get the character
////          ::std::cout << "tmp2: " << tmp2 << " ";
//          tmp ^= tmp2;        // zero the character in tmp
////          ::std::cout << "tmp0: " << tmp << " ";
//          tmp2 = static_cast<WORD_TYPE>(ALPHABET::TO_COMPLEMENT[tmp2 >> j]) << j;  // get the complement of the character
////          ::std::cout << "tmp2shift " << (tmp2 >> j) << " tmp2c: " << tmp2 << " ";
//
//          tmp ^= tmp2;  // replace the char
////          ::std::cout << "tmpc: " << tmp << ::std::endl;;
//        }
//
//        comp.data[i] = tmp;  // replace the word
//      }
//
//      // now reverse the whole thing sequentially.
//      bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SEQ>(result.getData(), comp.getData(), nWords);
//      result.right_shift_bits(nWords * sizeof(WORD_TYPE) * 8 - nBits);  // shift by remainder/padding.
//
//
//
////      // NOTE: left shift, pop, complement, then push + right shift is O(k * nWords), essentially O(k^2).
////      Kmer tmp_copy = src;
////      for (unsigned int i = 0; i < size; ++i)
////      {
////        result.do_left_shift(bitsPerChar);   // shifting the whole thing, inefficient but correct,
////                                            // especially for char that cross word boundaries.
////
////        // TODO: make this more efficient by operating on individual words.
////        tmp = ALPHABET::TO_COMPLEMENT[tmp_copy.data[0] & getLeastSignificantBitsMask<WORD_TYPE>(bitsPerChar)];
////
////        // copy `bitsperChar` least significant bits
////        copyBitsFixed<WORD_TYPE, bitsPerChar>(result.data[0], tmp);
////        tmp_copy.do_right_shift(bitsPerChar);
////      }
//
//      // result already was 0 to begin with, so no need to sanitize
//
//      return result;
//
//    }

    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA6>::value ||
                                    ::std::is_same<A, RNA6>::value ||
                                    ::std::is_same<A, DNA16>::value) &&
                                     ((sizeof(WORD_TYPE) * 8) % bitsPerChar == 0) &&
                                     ((sizeof(WORD_TYPE) * 8) < bitsPerChar), int>::type = 0>
    KMER_INLINE Kmer do_reverse_complement(Kmer const & src) const
    {

      /* Linear (inefficient) reverse:  because we don't have a specialization that supports TO_COMPLEMENT.  do word by word..*/

      // get temporary copy of this
      Kmer result(false);

      // get lower most bits from the temp copy and push them into the lower bits
      // of this
      WORD_TYPE tmp, tmp2, tmp3;
//      ::std::cout << ::std::hex;
      unsigned int charsPerWord = (sizeof(WORD_TYPE) * 8) / bitsPerChar;
      constexpr WORD_TYPE mask = getLeastSignificantBitsMask<WORD_TYPE>(bitsPerChar);

      // complement all, then reverse.
      for (unsigned int i = 0; i < nWords; ++i)
      {
        tmp = src.data[i];
        tmp3 = 0;

        for (unsigned int j = 0; j < charsPerWord; ++j) {
          tmp3 <<= bitsPerChar;
          tmp2 = tmp & mask;  // get the character
          tmp >>= bitsPerChar;
          tmp2 = static_cast<WORD_TYPE>(ALPHABET::to_complement(tmp2));  // get the complement of the character

          tmp3 |= tmp2;  // replace the char
        }

        result.data[nWords - 1 - i] = tmp3;  // replace the word
      }

      // now reverse the whole thing sequentially.
//      bliss::utils::bit_ops::reverse<bitsPerChar, ::bliss::utils::bit_ops::BIT_REV_SEQ>(result.getData(), comp.getData(), nWords);
      result.right_shift_bits(nWords * sizeof(WORD_TYPE) * 8 - nBits);  // shift by remainder/padding.


//      // NOTE: left shift, pop, complement, then push + right shift is O(k * nWords), essentially O(k^2).
//      Kmer tmp_copy = src;
//      for (unsigned int i = 0; i < size; ++i)
//      {
//        result.do_left_shift(bitsPerChar);   // shifting the whole thing, inefficient but correct,
//                                            // especially for char that cross word boundaries.
//
//        // TODO: make this more efficient by operating on individual words.
//        tmp = ALPHABET::TO_COMPLEMENT[tmp_copy.data[0] & getLeastSignificantBitsMask<WORD_TYPE>(bitsPerChar)];
//
//        // copy `bitsperChar` least significant bits
//        copyBitsFixed<WORD_TYPE, bitsPerChar>(result.data[0], tmp);
//        tmp_copy.do_right_shift(bitsPerChar);
//      }

      // result already was 0 to begin with, so no need to sanitize

      return result;

    }

    template <typename W = WORD_TYPE,
        typename ::std::enable_if<((sizeof(W) * 8) == bitsPerChar), int>::type = 0>
    KMER_INLINE Kmer do_reverse_complement(Kmer const & src) const
    {

      /* Linear (inefficient) reverse:  because we don't have a specialization that supports TO_COMPLEMENT.  do word by word..*/

      // get temporary copy of this
      Kmer result(false);

      // complement all, then reverse.
      for (unsigned int i = 0; i < nWords; ++i)
      {
        result.data[nWords - 1 - i] = static_cast<WORD_TYPE>(ALPHABET::to_complement(src.data[i]));  // replace the word
      }

      // now reverse the whole thing sequentially.
      result.right_shift_bits(nWords * sizeof(WORD_TYPE) * 8 - nBits);  // shift by remainder/padding.

      // result already was 0 to begin with, so no need to sanitize

      return result;

    }


  public:
    /**
     *
     * @brief   Create a new k-mer from the given sequence.
     *
     * The iterators base type has to be the same as the k-mers base type.
     * This function merely initializes the internal data with the values
     * from the iterator. There is no processing done on these values, which
     * especially means that there is no padding removed.
     *
     * The use of this function is mostly for testing purposes (creating the
     * expected reference (true) k-mers).
     *
     * @tparam InputIterator  An input iterator which is at least a forward
     *                        iterator.
     * @param begin           An interator pointing to the beginning of the
     *                        sequence to be used for initializing the k-mer
     *                        data structure.
     */
    template<typename InputIterator>
    explicit Kmer(InputIterator begin)
    {
      // assert that the iterator's value type is the same as this k-mers base
      // type
      static_assert(std::is_same<
          typename std::iterator_traits<InputIterator>::value_type,
          WORD_TYPE>::value,
          "Input iterator must have same value type as the Kmer storage");
  
      // copy all the data into this kmer
      WORD_TYPE* out = data;
      for (unsigned int i = 0; i < nWords; ++i)
      {
//        std::cout << " 0x" << std::hex << *begin;
        *(out++) = *(begin++);
      }
//      std::cout << std::endl;

      // set unused bits to zero to make this a valid kmer
      do_sanitize();
    }

  };


  /**
   * @brief << operator to write out DataBlock object's actual data.
   * @tparam Iterator   Source data iterator type.
   * @tparam Range      Range data type
   * @tparam Container  container type for buffer.  defaults to std::vector.
   * @param[in/out] ost   output stream to which the content is directed.
   * @param[in]     db    BufferedDataBlock object to write out
   * @return              output stream object
   */
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE=WordType>
  std::ostream& operator<<(std::ostream& ost, const Kmer<KMER_SIZE, ALPHABET, WORD_TYPE> & kmer)
  {
    ost << "Kmer=" << kmer.toString() << " (ASCII: " << bliss::utils::KmerUtils::toASCIIString(kmer) << ")";
  
    return ost;
  }
  

  } // namespace common
} // namespace bliss



#endif // BLISS_COMMON_KMER_H
