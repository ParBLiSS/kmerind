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
#include <stack>
#include <cstring>  // memset, memcpy


// own includes
#include "common/base_types.hpp"
#include "common/alphabets.hpp"
#include "common/alphabet_traits.hpp"
#include "common/bit_ops.hpp"
#include "common/padding.hpp"
#include "utils/kmer_utils.hpp"

// #include <xmmintrin.h>  // sse2
#if defined(__SSSE3__)
#include <tmmintrin.h>  // ssse3  includes sse2
#endif

#define INLINE inline

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

    using MACH_WORD_TYPE = size_t;

public:
    static constexpr int stride = sizeof(MACH_WORD_TYPE) / sizeof(WORD_TYPE);
    static constexpr int iters = nWords / stride;
    static constexpr int rem = nWords % stride;

    static constexpr int simd_stride = 16 / sizeof(WORD_TYPE);
    static constexpr int simd_iters = nWords / simd_stride;
    static constexpr int simd_rem = nWords % simd_stride;


  protected:
    /// The actual storage of the k-mer
    WORD_TYPE data[nWords]; //  not compatible with mxx datatype stuff - causes size of std::pair to increase greatly because alignment of elements is based on most resitrictive element.
                            // __attribute__((aligned(16)));
    						// actually, worse than that.  built-in types are having problems as well - alignment of primitive types is not reflected in sizeof or alignof keywords.
                            // TODO: can use MPI_Type_create_resized() to fix this problem for tuple and struct, but for primitive types its's not so easy.
  

    static constexpr uint8_t simd_mask_b[16] __attribute__((aligned(16))) = {0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F};
    static constexpr __m128i const *simd_mask = reinterpret_cast<const __m128i*>(simd_mask_b);


    // mask for reversing the order of items
    static constexpr uint8_t simd_rev_mask_b[16] __attribute__((aligned(16))) = {0x0F,0x0E,0x0D,0x0C,0x0B,0x0A,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,0x00};
    static constexpr __m128i const *simd_rev_mask = reinterpret_cast<const __m128i*>(simd_rev_mask_b);

    // lookup table for reversing bits in a byte in groups of 1
    static constexpr uint8_t simd_lut1_lo_b[16] __attribute__((aligned(16))) = {0x00,0x08,0x04,0x0c,0x02,0x0a,0x06,0x0e,0x01,0x09,0x05,0x0d,0x03,0x0b,0x07,0x0f};
    static constexpr __m128i const *simd_lut1_lo = reinterpret_cast<const __m128i*>(simd_lut1_lo_b);
    static constexpr uint8_t simd_lut1_hi_b[16] __attribute__((aligned(16))) = {0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0};
    static constexpr __m128i const *simd_lut1_hi = reinterpret_cast<const __m128i*>(simd_lut1_hi_b);
    // lookup table for reversing bits in a byte in groups of 2
    static constexpr uint8_t simd_lut2_lo_b[16] __attribute__((aligned(16))) = {0x00,0x04,0x08,0x0c,0x01,0x05,0x09,0x0d,0x02,0x06,0x0a,0x0e,0x03,0x07,0x0b,0x0f};
    static constexpr __m128i const *simd_lut2_lo = reinterpret_cast<const __m128i*>(simd_lut2_lo_b);
    static constexpr uint8_t simd_lut2_hi_b[16] __attribute__((aligned(16))) = {0x00,0x40,0x80,0xc0,0x10,0x50,0x90,0xd0,0x20,0x60,0xa0,0xe0,0x30,0x70,0xb0,0xf0};
    static constexpr __m128i const *simd_lut2_hi = reinterpret_cast<const __m128i*>(simd_lut2_hi_b);
    // lookup table for reversing bits in a byte in groups of 4 - no need, just use the simd_mask.

    // compare to self sets  register to all 1's.
    static constexpr uint8_t simd_true_mask_b[32] __attribute__((aligned(16))) = {0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
                                                                                  0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
    static constexpr __m128i const *simd_true_mask = reinterpret_cast<const __m128i*>(simd_true_mask_b);


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
    INLINE bool operator==(const Kmer& rhs) const
    {
      return memcmp(data, rhs.data, nWords * sizeof(WORD_TYPE)) == 0;
    }
  
    /**
     * @brief Compares this k-mer with the given k-mer for in-equality.
     *
     * Returns whether the data of the k-mers is different.
     *
     * @returns   `false` if the data of both k-mers is identical, `true`
     *            otherwise.
     */
    INLINE bool operator!=(const Kmer& rhs) const
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
    INLINE bool operator<(const Kmer& rhs) const
    {
      auto first = this->data;
      auto second = rhs.data;

      for (int i = 0; i < nWords; ++i, ++first, ++second) {
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
    INLINE bool operator<=(const Kmer& rhs) const
    {
      auto first = this->data;
      auto second = rhs.data;

      for (int i = 0; i < nWords; ++i, ++first, ++second) {
        if (*first != *second) return (*first < *second);  // if equal, keep comparing. else decide.
      }
      return true;  // all equal
//      std::pair<const WORD_TYPE*, const WORD_TYPE*> unequal = std::mismatch(this->data, this->data + nWords, rhs.data);
//      if (unequal.first == this->data + nWords)
//      {
//        // all elements are equal
//        return true;
//      }
//      else
//      {
//        // the comparison of the first unequal element will determine the result
//        return *(unequal.first) < *(unequal.second);
//      }
//      return memcmp(data, rhs.data, nWords * sizeof(WORD_TYPE)) <= 0;
    }
  
    /* symmetric comparison operators */
  
    /**
     * @brief Returns whether this k-mer compares greater or equal than the given
     *        k-mer.
     *
     * @returns `True` if this k-mer compares greater than or equal to the given
     *          k-mer, `False` otherwise.
     */
    INLINE bool operator>=(const Kmer& rhs) const
    {
      return ! (this->operator<(rhs));
    }
  
    /**
     * @brief Returns whether this k-mer compares greater than the given k-mer.
     *
     * @returns `True` if this k-mer compares greater than the given k-mer,
     *          `False` otherwise.
     */
    INLINE bool operator>(const Kmer& rhs) const
    {
      return ! (this->operator<=(rhs));
    }
  
    /* bit operators */
  
    /**
     * @brief XOR
     */
    INLINE Kmer& operator^=(const Kmer& rhs)
    {
      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_xor<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief XOR
     */
    INLINE Kmer operator^(const Kmer& rhs) const
    {
      Kmer result = *this;
      result ^= rhs;
      return result;
    }
  
    /**
     * @brief AND
     */
    INLINE Kmer& operator&=(const Kmer& rhs)
    {
      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_and<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief AND
     */
    INLINE Kmer operator&(const Kmer& rhs) const
    {
      Kmer result = *this;
      result &= rhs;
      return result;
    }
  
    /**
     * @brief OR
     */
    INLINE Kmer& operator|=(const Kmer& rhs)
    {
      std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_or<WORD_TYPE>());
      return *this;
    }
    /**
     * @brief OR
     */
    INLINE Kmer operator|(const Kmer& rhs) const
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
    INLINE Kmer& operator<<=(const std::size_t shift_by)
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
    INLINE Kmer operator<<(const std::size_t shift_by)
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
    INLINE Kmer& operator>>=(const std::size_t shift_by)
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
    INLINE Kmer operator>>(const std::size_t shift_by)
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
#if defined(__SSSE3__)
      return ::std::move(do_reverse_simd<ALPHABET>(*this));
#else
      return ::std::move(do_reverse_swar<ALPHABET>(*this));
#endif
    }

    /**
     * @brief Returns a reverse complement of a k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_complement() const
    {
#if defined(__SSSE3__)
      return ::std::move(do_reverse_complement_simd<ALPHABET>(*this));
#else
      return ::std::move(do_reverse_complement_swar<ALPHABET>(*this));
#endif
    }
  
  
    /**
     * @brief Returns a reversed k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_serial() const
    {
      return ::std::move(do_reverse_serial(*this));
    }

    /**
     * @brief Returns a reverse complement of a k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_complement_serial() const
    {
      return ::std::move(do_reverse_complement_serial(*this));
    }

#if defined(__SSSE3__)
    /**
     * @brief Returns a reversed k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_simd() const
    {
      return ::std::move(do_reverse_simd<ALPHABET>(*this));
    }

    /**
     * @brief Returns a reverse complement of a k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_complement_simd() const
    {
      return ::std::move(do_reverse_complement_simd<ALPHABET>(*this));
    }
#endif
    /**
     * @brief Returns a reversed k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_swar() const
    {
      return ::std::move(do_reverse_swar<ALPHABET>(*this));
    }

    /**
     * @brief Returns a reverse complement of a k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_complement_swar() const
    {
      return ::std::move(do_reverse_complement_swar<ALPHABET>(*this));
    }

    /**
     * @brief Returns a reversed k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_bswap() const
    {
      return ::std::move(do_reverse_bswap<ALPHABET>(*this));
    }

    /**
     * @brief Returns a reverse complement of a k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     * @returns   The reversed k-mer.
     */
    Kmer reverse_complement_bswap() const
    {
      return ::std::move(do_reverse_complement_bswap<ALPHABET>(*this));
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
  
    // for debug purposes
    /**
     * @brief Returns a string representation of this k-mer.
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
    INLINE void setBitsAtPos(WType w, int bitPos, int numBits) {
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
    INLINE void nextFromWordInternal(WType w, unsigned int shift)
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
    INLINE void nextReverseFromWordInternal(WType w, unsigned int shift)
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
    INLINE void do_sanitize()
    {
      // TODO use templated helper struct for <0> template specialization
      data[nWords-1] &= getLeastSignificantBitsMask<WORD_TYPE>(bitstream::invPadBits);
    }
  
    /**
     * @brief Sets all bits to zero.
     */
    INLINE void do_clear()
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
    INLINE void do_left_shift(size_t shift)
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
//      std::fill(data, data+word_shift, static_cast<WORD_TYPE>(0));
      memset(data, 0, word_shift * sizeof(WORD_TYPE));
    }
  
    /**
     * @brief Performs a right shift by `shift` bits on the k-mer data.
     *
     * @param shift   The number of bits to shift by.
     */
    // TODO template specialization for nWords = 1 (just use base type shift)
    // TODO implement more efficient version doing fixed left shift by BITS_PER_CHAR
    INLINE void do_right_shift(size_t shift)
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
//      std::fill(data + (nWords - word_shift), data + nWords, static_cast<WORD_TYPE>(0));
      memset(data + (nWords - word_shift), 0, word_shift * sizeof(WORD_TYPE));
    }
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    INLINE Kmer do_reverse_serial(Kmer const & src) const
    {
      // TODO implement logarithmic version (logarithmic in number of bits in a word, linear in number of words).  non power of 2 bitsPerChar is tricky because of byte boundaries.
  
      /* Linear (unefficient) reverse: */
  
      // get temporary copy of this
      Kmer tmp_copy = src;
      Kmer result;
  
      // get lower most bits from the temp copy and push them into the lower bits
      // of this
      for (unsigned int i = 0; i < size; ++i)
      {
        result.do_left_shift(bitsPerChar);   // shifting the whole thing, inefficient but correct,
                                            // especially for char that cross word boundaries.
        // copy `bitsperChar` least significant bits
        copyBitsFixed<WORD_TYPE, bitsPerChar>(result.data[0], tmp_copy.data[0]);
        tmp_copy.do_right_shift(bitsPerChar);
      }
  
      // result already was 0 to begin with, so no need to sanitize
//      // set ununsed bits to 0
//      this->do_sanitize();

      return ::std::move(result);
    }
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    INLINE Kmer do_reverse_complement_serial(Kmer const & src) const
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
        result.do_left_shift(bitsPerChar);   // shifting the whole thing, inefficient but correct,
                                            // especially for char that cross word boundaries.

        tmp = ALPHABET::TO_COMPLEMENT[tmp_copy.data[0] & getLeastSignificantBitsMask<WORD_TYPE>(bitsPerChar)];

        // copy `bitsperChar` least significant bits
        copyBitsFixed<WORD_TYPE, bitsPerChar>(result.data[0], tmp);
        tmp_copy.do_right_shift(bitsPerChar);
      }

      // result already was 0 to begin with, so no need to sanitize
//      // set ununsed bits to 0
//      this->do_sanitize();

      return ::std::move(result);

    }

#if defined(__SSSE3__)
    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.  8 to 100 times faster.
    // this code is inspired by http://stackoverflow.com/questions/746171/best-algorithm-for-bit-reversal-from-msb-lsb-to-lsb-msb-in-c
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, DNA>::value ||
    ::std::is_same<A, RNA>::value, __m128i>::type word_reverse_simd(__m128i const & b) {
      //== first shuffle (reverse) the bytes
      __m128i r = _mm_shuffle_epi8(b, *simd_rev_mask);                                    // SSSE3

      //== get the lut indices.
      // lower 4 bits
      __m128i lo = _mm_and_si128(*simd_mask, r);                                          // SSE2
      // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
      __m128i hi = _mm_srli_epi16(_mm_andnot_si128(*simd_mask, r), 4);                    // SSE2

      //== now use shuffle.  hi and lo are now indices. and lut2_lo/hi are lookup tables.  remember that hi/lo are swapped.
      return _mm_or_si128(_mm_shuffle_epi8(*simd_lut2_hi, lo), _mm_shuffle_epi8(*simd_lut2_lo, hi));  // SSSE3, SSE2

    }
    /// do reverse complement.  8 to 100 times faster than serially reversing the bits.
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, DNA>::value ||
        ::std::is_same<A, RNA>::value, __m128i>::type word_reverse_complement_simd(__m128i const & b) {
      // DNA type, requires that alphabet be setup so that negation produces the complement.
      return _mm_xor_si128(word_reverse_simd<A>(b), *simd_true_mask);  // sse does not have negation.   // SSE2
    }

    // reverse via bit swapping in 4 bit increment.  this complement table is then used for lookup.
    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.
    /// 8 to 50 times faster than sequentially reverse the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, DNA16>::value, __m128i>::type word_reverse_simd(__m128i const & b) {
      //== first shuffle (reverse) the bytes
      __m128i r = _mm_shuffle_epi8(b, *simd_rev_mask);

      // lower 4 bits shift to upper
      __m128i lo = _mm_slli_epi32(_mm_and_si128(*simd_mask, r), 4);
      // upper 4 bits shift to lower.
      __m128i hi = _mm_srli_epi32(_mm_andnot_si128(*simd_mask, r), 4);

      // swap bits in groups of 4 bits, so after this point, just or.

      return _mm_or_si128(hi, lo);
    }
    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET>
    static INLINE typename ::std::enable_if<::std::is_same<A, DNA16>::value, __m128i>::type word_reverse_complement_simd(__m128i const & b) {
      //== first shuffle (reverse) the bytes
      __m128i r = _mm_shuffle_epi8(b, *simd_rev_mask);

      //== get the lut indices.
      // lower 4 bits
      __m128i lo = _mm_and_si128(*simd_mask, r);
      // upper 4 bits.  note that after mask, we shift by 4 bits int by int.  each int will have the high 4 bits set to 0 from the shift, but they were going to be 0 anyways.
      __m128i hi = _mm_srli_epi32(_mm_andnot_si128(*simd_mask, r), 4);

      //== now use shuffle to look up the reversed bytes.  hi and lo are now indices. and lut2_lo/hi are lookup tables. remember that lo and hi need to be swapped.
      return _mm_or_si128(_mm_shuffle_epi8(*simd_lut1_hi, lo), _mm_shuffle_epi8(*simd_lut1_lo, hi));

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
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    INLINE Kmer do_reverse_simd(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output

      // do the rem.  When rem is not 0, then it's very slow.  (like 40 to 50 times slower.)  (because of maskmoveu)
      if (simd_rem > 0) {
        __m128i mword = word_reverse_simd<A>(_mm_loadu_si128((__m128i*)(src.data + nWords - simd_stride)));     // SSE2
        // store back to result - uses a mask to select which bytes are stored.
        if (simd_iters > 0) {
          // has room, extra will be overwritten, so just use _mm_storeu.
            _mm_storeu_si128((__m128i*)(result.data), mword);                                                   // SSE2
        } else {
          // not enough room, so need to copy.  maskmoveu is very slow.  why?
          memcpy(result.data, &mword, simd_rem * sizeof(WORD_TYPE));
            //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));               // SSE2
        }
      }  // else no rem.


      // swap the word order.  also swap the packed chars in the word at the same time.
      if (simd_iters > 0) {
        auto in = src.data;
        auto out = result.data + nWords - simd_stride;

        for (int i = 0; i < simd_iters; ++i) {
          // load, reverse, then store.  assuming unaligned.
          _mm_storeu_si128((__m128i*)out, word_reverse_simd<A>(_mm_loadu_si128((__m128i*)in)));                 // SSE2
          in += simd_stride;
          out -= simd_stride;
        }
      }


      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return ::std::move(result);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    INLINE Kmer do_reverse_complement_simd(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output

      // do the rem.  When rem is not 0, then it's very slow.  (like 40 to 50 times slower.)  (because of maskmoveu)
      if (simd_rem > 0) {
        __m128i mword = word_reverse_complement_simd<A>(_mm_loadu_si128((__m128i*)(src.data + nWords - simd_stride)));   // SSE2
        // store back to result - uses a mask to select which bytes are stored.
        // do the remainder first because we want the reversed high bits to be overwritten.
        if (simd_iters > 0) {
          // has room, extra will be overwritten, so just use _mm_storeu.
            _mm_storeu_si128((__m128i*)(result.data), mword);                                                             // SSE2
        } else {
          // not enough room, so do maskmoveu.  maskmoveu is very slow.  why?
// don't need this part...
//          WORD_TYPE tmp[simd_stride] __attribute__((aligned(16)));
//          _mm_storeu_si128((__m128i*)tmp, mword);
          memcpy(result.data, &mword, simd_rem * sizeof(WORD_TYPE));

          //_mm_maskmoveu_si128(mword, *simd_store_mask, reinterpret_cast<char*>(result.data));                           // SSE2
        }
      }  // else no rem.


      // swap the word order.  also swap the packed chars in the word at the same time.
      if (simd_iters > 0) {
        auto in = src.data;
        auto out = result.data + nWords - simd_stride;

        for (int i = 0; i < simd_iters; ++i) {
          // load, reverse, then store.  assuming unaligned.
          _mm_storeu_si128((__m128i*)out, word_reverse_complement_simd<A>(_mm_loadu_si128((__m128i*)in)));                // SSE2
          in += simd_stride;
          out -= simd_stride;
        }
      }



      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return ::std::move(result);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA16>::value), int>::type = 0>
    INLINE Kmer do_reverse_simd(Kmer const & src) const
    {
       return ::std::move(do_reverse_serial(src));
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA16>::value), int>::type = 0>
    INLINE Kmer do_reverse_complement_simd(Kmer const & src) const
    {
      return ::std::move(do_reverse_complement_serial(src));
    }


#endif

    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.  8 to 100 times faster.
    template <typename A = ALPHABET>
    INLINE typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                     ::std::is_same<A, RNA>::value, WORD_TYPE>::type word_reverse(WORD_TYPE const & b) const {
      WORD_TYPE v = b;
      WORD_TYPE mask = ~0;
      mask ^= (mask << (sizeof(WORD_TYPE) * 4));

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      switch (sizeof(WORD_TYPE) * 8)
      {
        // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
        case 64 :  v =  (v >> 32)         |  (v         << 32);  mask ^= (mask << 16);
        case 32 :  v = ((v >> 16) & mask) | ((v & mask) << 16);  mask ^= (mask << 8);
        case 16 :  v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);
        case 8  :  v = ((v >>  4) & mask) | ((v & mask) <<  4);  mask ^= (mask << 2);
                   v = ((v >>  2) & mask) | ((v & mask) <<  2);
        default :  break;
      }

      return v;
    }
    /// do reverse complement.  8 to 100 times faster than serially reversing the bits.
    template <typename A = ALPHABET>
        INLINE typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                         ::std::is_same<A, RNA>::value, WORD_TYPE>::type word_reverse_complement(WORD_TYPE const & b) const {
      // DNA type, requires that alphabet be setup so that negation produces the complement.
      return ~(word_reverse<A>(b));
    }

    // reverse via bit swapping in 4 bit increment.  this complement table is then used for lookup.
    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.
    /// 8 to 50 times faster than sequentially reverse the bits.
    template <typename A = ALPHABET>
        INLINE typename ::std::enable_if<::std::is_same<A, DNA16>::value, WORD_TYPE>::type  word_reverse(WORD_TYPE const & b) const {
      WORD_TYPE v = b;
      WORD_TYPE mask = ~0;
      mask ^= (mask << (sizeof(WORD_TYPE) * 4));

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.
      switch (sizeof(WORD_TYPE) * 8)
      {
        // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
        case 64 :  v =  (v >> 32)         |  (v         << 32);  mask ^= (mask << 16); // can use BSWAP_64 here
        case 32 :  v = ((v >> 16) & mask) | ((v & mask) << 16);  mask ^= (mask << 8);  // can use BSWAP here.
        case 16 :  v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);  // can use XCHG here
        case 8  :  v = ((v >>  4) & mask) | ((v & mask) <<  4);
        default :  break;
      }

      return v;
    }
    /// do reverse complement.  8 to 50 times faster than serially reversing the bits.
    template <typename A = ALPHABET>
            INLINE typename ::std::enable_if<::std::is_same<A, DNA16>::value, WORD_TYPE>::type word_reverse_complement(WORD_TYPE const & b) const {
      // DNA type:
      WORD_TYPE v = b;
      WORD_TYPE mask = ~0;
      mask ^= (mask << (sizeof(WORD_TYPE) * 4));

      // essentially a Duff's Device here - unrolled loop with constexpr to allow compiler to optimize here.

      // REQUIRES that alphabet where complement is bit reversal of the original character.
      // since for DNA16, complement is bit reversal in a nibble, we do 2 extra iterations.
      switch (sizeof(WORD_TYPE) * 8)
      {
        // TODO: do byte level shuffling in 1 instruction. - would need to take care of WORD_TYPE vs machine word size mismatch.
        case 64 :  v =  (v >> 32)         |  (v         << 32);  mask ^= (mask << 16);
        case 32 :  v = ((v >> 16) & mask) | ((v & mask) << 16);  mask ^= (mask << 8);
        case 16 :  v = ((v >>  8) & mask) | ((v & mask) <<  8);  mask ^= (mask << 4);
        case 8  :  v = ((v >>  4) & mask) | ((v & mask) <<  4);  mask ^= (mask << 2);  // reverse nibbles within a byte
                   v = ((v >>  2) & mask) | ((v & mask) <<  2);  mask ^= (mask << 1);  // last 2 steps create the complement
                   v = ((v >>  1) & mask) | ((v & mask) <<  1);
        default :  break;
      }

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
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    INLINE Kmer do_reverse_swar(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output

      // swap the word order.  also swap the packed chars in the word at the same time.
      for (int i = 0, j = nWords - 1, max = nWords; i < max; ++i, --j)
        result.data[j] = word_reverse<A>(src.data[i]);

      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return ::std::move(result);
    }
  
    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    INLINE Kmer do_reverse_complement_swar(Kmer const & src) const
    {
      Kmer result;  // empty kmerfor output

      // swap the word order.  also do rev_comp for the chars in the word at the same time.
      for (int i = 0, j = nWords - 1, max = nWords; i < max; ++i, --j)
        result.data[j] = word_reverse_complement<A>(src.data[i]);

      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return ::std::move(result);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA16>::value), int>::type = 0>
    INLINE Kmer do_reverse_swar(Kmer const & src) const
    {
       return ::std::move(do_reverse_serial(src));
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA16>::value), int>::type = 0>
    INLINE Kmer do_reverse_complement_swar(Kmer const & src) const
    {
      return ::std::move(do_reverse_complement_serial(src));
    }



    // ================ experiement with bswap operations.
#include <bits/byteswap.h>

    static constexpr uint64_t mask0 = 0x0F0F0F0F0F0F0F0F;
    static constexpr uint64_t mask1 = 0x3333333333333333;
    static constexpr uint64_t mask2 = 0x5555555555555555;
//    static constexpr uint32_t mask_32[3] = { 0x0F0F0F0F, 0x33333333, 0x55555555 };
//    static constexpr uint16_t mask_16[3] = { 0x0F0F, 0x3333, 0x5555 };
//    static constexpr uint8_t  mask_8[3]  = { 0x0F, 0x33, 0x55 };

    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.  8 to 100 times faster.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
    INLINE typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                     ::std::is_same<A, RNA>::value, W>::type word_reverse_bswap(W const & b) const {
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
        INLINE typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                         ::std::is_same<A, RNA>::value, W>::type word_reverse_complement_bswap(W const & b) const {
      // DNA type, requires that alphabet be setup so that negation produces the complement.
      return ~(word_reverse_bswap<A, W>(b));
    }

    // reverse via bit swapping in 4 bit increment.  this complement table is then used for lookup.
    /// reverse packed characters in a word. implementation is compatible with alphabet size of 1, 2 or 4 bits.
    /// 8 to 50 times faster than sequentially reverse the bits.
    // TODO:  the use of bswap_xx makes it gnu c dependent.
    template <typename A = ALPHABET, typename W = WORD_TYPE>
        INLINE typename ::std::enable_if<::std::is_same<A, DNA16>::value, W>::type  word_reverse_bswap(W const & b) const {
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
        INLINE typename ::std::enable_if<::std::is_same<A, DNA16>::value, W>::type  word_reverse_complement_bswap(W const & b) const {
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
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    INLINE Kmer do_reverse_bswap(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output

      // unchunked reverse works.
            for (int i = 0, j = nWords - 1, max = nWords; i < max; ++i, --j)
              result.data[j] = word_reverse_bswap<A, WORD_TYPE>(src.data[i]);


/*      //== chunked revcomp does NOT always work, although slightly faster, especially for small words.
      if (iters > 0) {
    	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
              const MACH_WORD_TYPE *s = reinterpret_cast<const MACH_WORD_TYPE*>(src.data);
              WORD_TYPE *d = result.data + nWords - stride;

              // swap the word order.  also swap the packed chars in the word at the same time.
              // do the first part of source in multiple of MACH_WORD_TYPE.
              for (int i = 0; i < iters; ++i, d -= stride)
                *(reinterpret_cast<MACH_WORD_TYPE*>(d)) = word_reverse_bswap<A, MACH_WORD_TYPE>(s[i]);
            }
      //      printf("1st step:\t\t%s\n", result.toAlphabetString().c_str());

            // do the rem.
            if (rem > 1) {
          	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
              MACH_WORD_TYPE mword = word_reverse_bswap<A, MACH_WORD_TYPE>(*(reinterpret_cast<const MACH_WORD_TYPE*>(src.data + nWords - stride)));
              memcpy(result.data, &mword, rem * sizeof(WORD_TYPE));
            } else if (rem == 1) {
              result.data[0] = word_reverse_bswap<A, WORD_TYPE>(src.data[nWords - 1]);
            }  // else no rem.
      //      printf("2nd step:\t\t%s\n", result.toAlphabetString().c_str());
*/

            // shift if necessary      // ununsed bits will be set to 0 by shift
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return ::std::move(result);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<::std::is_same<A, DNA>::value ||
                                  ::std::is_same<A, RNA>::value ||
                                  ::std::is_same<A, DNA16>::value, int>::type = 0>
    INLINE Kmer do_reverse_complement_bswap(Kmer const & src) const
    {
      Kmer result;  // empty kmer for output


      // swap the word order.  also do rev_comp for the chars in the word at the same time.
      for (int i = 0, j = nWords - 1, max = nWords; i < max; ++i, --j)
        result.data[j] = word_reverse_complement_bswap<A, WORD_TYPE>(src.data[i]);

/*      //== chunked revcomp does NOT always work, although slightly faster, especially for small words.
      if (iters > 0) {
    	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
        const MACH_WORD_TYPE *s = reinterpret_cast<const MACH_WORD_TYPE*>(src.data);
        WORD_TYPE *d = result.data + nWords - stride;

        // swap the word order.  also swap the packed chars in the word at the same time.
        // do the first part of source in multiple of MACH_WORD_TYPE.
        for (int i = 0; i < iters; ++i, d -= stride)
          *(reinterpret_cast<MACH_WORD_TYPE*>(d)) = word_reverse_complement_bswap<A, MACH_WORD_TYPE>(s[i]);
      }
//      printf("1st step:\t\t%s\n", result.toAlphabetString().c_str());

      // do the rem.
      if (rem > 1) {
    	  // right here there is a warrning about 'dereferencing type-punned pointer will break strict-aliasing rules'
        MACH_WORD_TYPE mword = word_reverse_complement_bswap<A, MACH_WORD_TYPE>(*(reinterpret_cast<const MACH_WORD_TYPE*>(src.data + nWords - stride)));
        memcpy(result.data, &mword, rem * sizeof(WORD_TYPE));
      } else if (rem == 1) {
        result.data[0] = word_reverse_complement_bswap<A, WORD_TYPE>(src.data[nWords - 1]);
      }  // else no rem.
//      printf("2nd step:\t\t%s\n", result.toAlphabetString().c_str());
*/
      // shift if necessary      // ununsed bits will be set to 0 by shift
      result.do_right_shift(nWords * sizeof(WORD_TYPE) * 8 - nBits);

      return ::std::move(result);
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA16>::value), int>::type = 0>
    INLINE Kmer do_reverse_bswap(Kmer const & src) const
    {
       return ::std::move(do_reverse_serial(src));
    }

    /**
     * @brief Reverses this k-mer.
     *
     * Note that this does NOT reverse the bit pattern, but reverses
     * the sequence of `BITS_PER_CHAR` bits each.
     *
     */
    template <typename A = ALPHABET,
        typename ::std::enable_if<!(::std::is_same<A, DNA>::value ||
                                    ::std::is_same<A, RNA>::value ||
                                    ::std::is_same<A, DNA16>::value), int>::type = 0>
    INLINE Kmer do_reverse_complement_bswap(Kmer const & src) const
    {
      return ::std::move(do_reverse_complement_serial(src));
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
  
  
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_mask_b[16];
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_rev_mask_b[16];
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_lut1_lo_b[16];
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_lut1_hi_b[16];
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_lut2_lo_b[16];
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_lut2_hi_b[16];
  template<unsigned int KMER_SIZE, typename ALPHABET, typename WORD_TYPE>
  constexpr uint8_t Kmer<KMER_SIZE, ALPHABET, WORD_TYPE>::simd_true_mask_b[32];


  } // namespace common
} // namespace bliss



#endif // BLISS_COMMON_KMER_H
