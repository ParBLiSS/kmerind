/**
 * @file    kmer.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements the Kmer data type.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
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
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>  // std::pair


// own includes
#include <common/base_types.hpp>
#include <common/alphabets.hpp>
#include <common/alphabet_traits.hpp>
#include <common/bit_ops.hpp>
#include <common/padding.hpp>
#include <utils/kmer_utils.hpp>


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
  
  protected:
    /// The actual storage of the k-mer
    WORD_TYPE data[nWords];
  
  public:
  
    /**
     * @brief   The default constructor, creates an uninitialized k-mer.
     *
     * The k-mer's data is not initalized, but still represents a valid k-mer,
     * as the unused bits are set to 0.
     */
    Kmer()
    {
      // make this a valid Kmer, other bits are left uninitialized
      // do_sanitize();
  
      // set all bits to zero (no more 'may be used uninitialized' warning)
      do_clear();
    }
  
    /// copy constructor
    Kmer(const Kmer& other) : Kmer(other.data) {};
  
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
      for (int i = 0; i < KMER_SIZE; ++i, bitPos -= bitsPerChar)
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
      for (int i = 0; i < KMER_SIZE; ++i, bitPos += bitsPerChar)
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
      for (int i = 0; i < KMER_SIZE; ++i, bitPos -= bitsPerChar) {
  
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
      for (int i = 0; i < KMER_SIZE; ++i, bitPos += bitsPerChar) {
  
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
      //do_sanitize();
  
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
  //    do_sanitize();
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
    inline bool operator==(const Kmer& rhs) const
    {
      return std::equal(data, data+nWords, rhs.data);
    }
  
    /**
     * @brief Compares this k-mer with the given k-mer for in-equality.
     *
     * Returns whether the data of the k-mers is different.
     *
     * @returns   `false` if the data of both k-mers is identical, `true`
     *            otherwise.
     */
    inline bool operator!=(const Kmer& rhs) const
    {
      return !(this->operator==(rhs));
    }
  
    /* ordered comparison operators */
  
    /**
     * @brief Returns whether this k-mer compares smaller than the given k-mer.
     *
     * @returns `True` if this k-mer compares smaller than the given k-mer,
     *          `False` otherwise.
     */
    inline bool operator<(const Kmer& rhs) const
    {
      std::pair<const WORD_TYPE*, const WORD_TYPE*> unequal = std::mismatch(this->data, this->data + nWords, rhs.data);
      if (unequal.first == this->data + nWords)
      {
        // all elements are equal
        return false;
      }
      else
      {
        // the comparison of the first unequal element will determine the result
        return *(unequal.first) < *(unequal.second);
      }
    }
  
    /**
     * @brief Returns whether this k-mer compares smaller or equal than the given
     *        k-mer.
     *
     * @returns `True` if this k-mer compares smaller than or equal to the given
     *          k-mer, `False` otherwise.
     */
    inline bool operator<=(const Kmer& rhs) const
    {
      std::pair<const WORD_TYPE*, const WORD_TYPE*> unequal = std::mismatch(this->data, this->data + nWords, rhs.data);
      if (unequal.first == this->data + nWords)
      {
        // all elements are equal
        return true;
      }
      else
      {
        // the comparison of the first unequal element will determine the result
        return *(unequal.first) < *(unequal.second);
      }
    }
  
    /* symmetric comparison operators */
  
    /**
     * @brief Returns whether this k-mer compares greater or equal than the given
     *        k-mer.
     *
     * @returns `True` if this k-mer compares greater than or equal to the given
     *          k-mer, `False` otherwise.
     */
    inline bool operator>=(const Kmer& rhs) const
    {
      return ! (this->operator<(rhs));
    }
  
    /**
     * @brief Returns whether this k-mer compares greater than the given k-mer.
     *
     * @returns `True` if this k-mer compares greater than the given k-mer,
     *          `False` otherwise.
     */
    inline bool operator>(const Kmer& rhs) const
    {
      return ! (this->operator<=(rhs));
    }
  
    /* bit operators */
  
    /**
     * @brief XOR
     */
    inline Kmer& operator^=(const Kmer& rhs)
    {
      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_xor<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief XOR
     */
    inline Kmer operator^(const Kmer& rhs) const
    {
      Kmer result = *this;
      result ^= rhs;
      return result;
    }
  
    /**
     * @brief AND
     */
    inline Kmer& operator&=(const Kmer& rhs)
    {
      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_and<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief AND
     */
    inline Kmer operator&(const Kmer& rhs) const
    {
      Kmer result = *this;
      result &= rhs;
      return result;
    }
  
    /**
     * @brief OR
     */
    inline Kmer& operator|=(const Kmer& rhs)
    {
      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_or<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief OR
     */
    inline Kmer operator|(const Kmer& rhs) const
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
    inline Kmer& operator<<=(const std::size_t shift_by)
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
    inline Kmer operator<<(const std::size_t shift_by)
    {
      Kmer result = *this;
      result <<= shift_by;
      return result;
    }
  
    /**
     * @brief Shifts the k-mer right by the given number of CHARACTERS.
     *
     * @note  This shifts by the number of characters (which is a larger shift
     *        then bitwise).
     */
    inline Kmer& operator>>=(const std::size_t shift_by)
    {
      // shift the given number of **characters** !
      this->do_right_shift(shift_by * bitsPerChar);
      this->do_sanitize();
      return *this;
    }
  
  
    /**
     * @brief Shifts the k-mer right by the given number of CHARACTERS.
     *
     * @note  This shifts by the number of characters (which is a larger shift
     *        then bitwise).
     */
    inline Kmer operator>>(const std::size_t shift_by)
    {
      Kmer result = *this;
      result >>= shift_by;
      return result;
    }
  
  
    /**
     * @brief Returns a reversed k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reversed_kmer() const
    {
      // create a copy of this
      Kmer rev(*this);
      // reverse it
      rev.do_reverse();
      // return the result
      return rev;
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
      ss << "k-mer of size " << size << ": [";
      for (unsigned int i = 0; i < nWords; ++i)
      {
        ss << "0x" << std::hex << data[i] << " ";
      }
      ss << "]";
      return ss.str();
    }
  
  
    std::string toAlphabetString() const
    {
      std::stringstream ss;
  
      WordType w = 0;
      int bit_pos = 0;
      int word_pos = 0;
      int bit_pos_in_word = 0;
      int offset = 0;
      for (int i = KMER_SIZE - 1; i >= 0; --i)
      {
        // start pos
        bit_pos = i * bitsPerChar;
        word_pos = bit_pos / (sizeof(WordType) * 8);
        bit_pos_in_word = bit_pos % (sizeof(WordType) * 8);
  
        w = (data[word_pos] >> bit_pos_in_word);
  
        if (bit_pos_in_word + bitsPerChar > (sizeof(WordType) * 8)) {  // == means that there were just enough bits in the previous word.
          offset = (sizeof(WordType) * 8) - bit_pos_in_word;
          ++word_pos;
  
          // now need to get the next part.
          if (word_pos < nWords) w = w | (data[word_pos] << offset);
  
        }
        w = w & getLeastSignificantBitsMask<WORD_TYPE>(bitsPerChar);
  
  
        ss << w << " ";
  
      }
  
      return ss.str();
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
        for (int i = nwords - 1; i >= 0; --i) {
  
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
    inline void setBitsAtPos(WType w, int bitPos, int numBits) {
      // get the value masked.
      WORD_TYPE charVal = static_cast<WORD_TYPE>(w) & getLeastSignificantBitsMask<WORD_TYPE>(numBits);
  
      // determine which word in kmer it needs to go to
      int wordId = bitPos / bitstream::bitsPerWord;
      int offsetInWord = bitPos % bitstream::bitsPerWord;  // offset is where the LSB of the char will sit, in bit coordinate.
  
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
    inline void nextFromWordInternal(WType w, unsigned int shift)
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
    inline void nextReverseFromWordInternal(WType w, unsigned int shift)
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
    inline void do_sanitize()
    {
      // TODO use templated helper struct for <0> template specialization
      data[nWords-1] &= getLeastSignificantBitsMask<WORD_TYPE>(bitstream::invPadBits);
    }
  
    /**
     * @brief Sets all bits to zero.
     */
    inline void do_clear()
    {
      std::fill(data, data + nWords, static_cast<WORD_TYPE>(0));
    }
  
    /**
     * @brief Performs a left shift by `shift` bits on the k-mer data.
     *
     * @param shift   The number of bits to shift by.
     */
    // TODO template specialization for nWords = 1 (just use base type shift)
    // TODO implement more efficient version doing fixed left shift by BITS_PER_CHAR
    inline void do_left_shift(size_t shift)
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
        for (size_t i = nWords - 1; i >= word_shift; --i)
        {
          data[i] = data[i - word_shift];
        }
      }
      else
      {
        const size_t inv_offset = sizeof(WORD_TYPE)*8 - offset;
        for (size_t i = nWords - 1; i > word_shift; --i)
        {
          data[i] = ((data[i - word_shift] << offset) | (data[i - word_shift - 1] >> inv_offset));
        }
        data[word_shift] = data[0] << offset;
      }
  
      // set all others to 0
      std::fill(data, data+word_shift, static_cast<WORD_TYPE>(0));
    }
  
    /**
     * @brief Performs a right shift by `shift` bits on the k-mer data.
     *
     * @param shift   The number of bits to shift by.
     */
    // TODO template specialization for nWords = 1 (just use base type shift)
    // TODO implement more efficient version doing fixed left shift by BITS_PER_CHAR
    inline void do_right_shift(size_t shift)
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
        for (size_t i = 0; i < nWords - word_shift - 1; ++i)
        {
          data[i] = ((data[i + word_shift] >> offset) | (data[i + word_shift + 1] << inv_offset));
        }
        data[nWords - word_shift - 1] = data[nWords-1] >> offset;
      }
      // set all others to 0
      std::fill(data + (nWords - word_shift), data + nWords, static_cast<WORD_TYPE>(0));
    }
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    inline void do_reverse()
    {
      // TODO implement logarithmic version (logarithmic in number of bits)
  
      /* Linear (unefficient) reverse: */
  
      // get temporary copy of this
      Kmer tmp_copy = *this;
  
      // get lower most bits from the temp copy and push them into the lower bits
      // of this
      for (unsigned int i = 0; i < size; ++i)
      {
        this->do_left_shift(bitsPerChar);   // shifting the whole thing, inefficient but correct,
                                            // especially for char that cross word boundaries.
        // copy `bitsperChar` least significant bits
        copyBitsFixed<WORD_TYPE, bitsPerChar>(this->data[0], tmp_copy.data[0]);
        tmp_copy.do_right_shift(bitsPerChar);
      }
  
      // set ununsed bits to 0
      this->do_sanitize();
    }
  
  
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
        *(out++) = *(begin++);
      }
  
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
