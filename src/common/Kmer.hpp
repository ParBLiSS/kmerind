/**
 * @file    Kmer.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Implements the Kmer data type.
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
#include <string>
#include <sstream>
#include <algorithm>
#include <utility>

// own includes
#include <common/base_types.hpp>
#include <common/bit_ops.hpp>
#include <common/padding.hpp>

namespace bliss
{

/**
 * @brief   Implements a general templated k-mer class.
 *
 * The k-mer size (number of characters) is a template parameter, thus
 * the k-mer size for this k-mer needs to be fixed at compile time. The given
 * k-mer size and the `BITS_PER_CHAR` determine the size of the underlying
 * data array at compile time.
 *
 * @todo: implement and refer to the dynamic k-mer (which will be more
 *        inefficient)
 *
 * @tparam  KMER_SIZE   The size (number of characters) of the k-mer.
 * @tparam  BITS_PER_CHAR   The number of bits per character in the k-mer.
 *                          E.g. for DNA this would be 2.
 * @tparam  word_type       The unsigned integer type to be used as the
 *                          storage type. This can be of any size from
 *                          uint8_t up to uint64_t. (default is uint64_t)
 */
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

  /*
   * last character offsets (and whether or not it is split accord storage
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

  /**
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
  Kmer(InputIterator begin)
  {
    // assert that the iterator's value type is the same as this k-mers base
    // type
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
   *  - hash function (!?)
   */

  /**
   * FIXME: update documentation for `size => (size-1)`
   * @brief   Fills this k-mer from a packed and padded input sequence.
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
  template <typename InputIterator>
  // TODO: add option for bit offset in input sequence?
  unsigned int fillFromPaddedStream(InputIterator& begin, bool stop_on_last = false)
  {
    // remove padding and copy to own data structure
    typedef typename std::iterator_traits<InputIterator>::value_type input_word_type;
    const unsigned int paddingBits = PaddingTraits<input_word_type, bitsPerChar>::padding_bits;
    removePadding(begin, data, size*bitsPerChar, paddingBits);

    // set unused bits to 0
    do_sanitize();

    // the bit offset in the input sequence
    // TODO: do this inside the removePadding function!?
    unsigned int total_bits = stop_on_last ? (size-1)*bitsPerChar : size*bitsPerChar;
    unsigned int offset = total_bits % PaddingTraits<input_word_type, bitsPerChar>::data_bits;
    std::advance(begin, total_bits / PaddingTraits<input_word_type, bitsPerChar>::data_bits);

    // return the offset
    return offset;
  }


  template <typename InputIterator>
  void fillFromChars(InputIterator& begin, bool stop_on_last = false)
  {
    // value type of given iterator
    typedef typename std::iterator_traits<InputIterator>::value_type char_type;

    // clear k-mer
    do_clear();

    // add to lsb iteratively and reverse at the end
    for (unsigned int i = 0; i < size;)
    {
      // get next character as word_type and mask out bits that are not needed
      char_type c = *begin;
      word_type w = static_cast<word_type>(c);
      w &= getBitMask<word_type>(bitsPerChar);

      // left shift k-mer
      // TODO: replace by single shift operation
      do_left_shift(bitsPerChar);

      // add character to least significant end (requires least shifting)
      *data |= w;

      // iterate the input iterator one more, but stop on last
      // iteration if that option is set
      if (!stop_on_last || (stop_on_last && ++i!=size)) ++begin;
    }

    // reverse the k-mer (needed since we added the newest characters to the
    // least significant end rather than the most significant end)
    reversed_kmer();
  }

  /**
   * @brief   Generates the next k-mer from the given sequence using a sliding
   *          window approach.
   *
   * Given the packed and padded input sequence via the `begin` iterator
   * and the current bit offset, this function will read the next
   * `BITS_PER_CHAR` bits, right shift the current k-mer value by that number
   * of bits and puts the newly read bits into the `BITS_PER_CHAR` most
   * significant bits of the k-mer value. Both parameters are passed by
   * reference and internally updated to the next position to read from.
   * Therefore, the user of this function just needs to keep calling this
   * function without updating the parameters to generate all k-mers for a
   * given packed and padded input sequence.
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
  void nextFromPaddedStream(InputIterator& begin, offset_t& offset)
  {
    typedef typename std::iterator_traits<InputIterator>::value_type input_type;
    // shift the kmer by the size of one character
    // TODO: replace this by a call that does exactly bitsPerChar right shift
    // (better compiler optimization)
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
    if (offset >= PaddingTraits<input_type, bitsPerChar>::data_bits)
    {
      ++begin;
      offset -= PaddingTraits<input_type, bitsPerChar>::data_bits;
    }
  }

  /**
   * @brief Creates the next k-mer by shifting and adding the given character.
   *
   * Creates the next k-mer by the sliding window. This first right shifts the
   * current k-mer and then adds the given character to the most significant
   * bits of this k-mer.
   *
   * Note: this is changing this k-mer in-place.
   *
   * @param c The character to be added to the k-mer.
   */
  void nextFromChar(char c)
  {
    // shift the kmer by the size of one character
    // TODO: replace this by a call that does exactly bitsPerChar right shift
    // (better compiler optimization)
    do_right_shift(bitsPerChar);

    // cast the char into our word size and then AND it with a bitmask
    word_type new_char = static_cast<word_type>(c);
    new_char &= getBitMask<word_type>(bitsPerChar);

    // shift it to the right position and OR it into our data
    int shift_by = lastCharWordOffset;
    data[nWords - 1 - lastCharIsSplit] |= static_cast<word_type>(new_char << shift_by);

    // the last character might be split between storage words
    // NOTE: this branch is removed by the compiler (result known at compile time)
    if (lastCharIsSplit)
    {
      data[nWords - 1] = static_cast<word_type>(new_char >> leftSplitSize);
    }

    // clean up
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
    std::pair<const word_type*, const word_type*> unequal = std::mismatch(this->data, this->data + nWords, rhs.data);
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
    std::pair<const word_type*, const word_type*> unequal = std::mismatch(this->data, this->data + nWords, rhs.data);
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
    std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_xor<word_type>());
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
    std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_and<word_type>());
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
    std::transform(this->data, this->data + nWords, rhs.data, this->data, std::bit_or<word_type>());
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
    return *this;
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
    return *this;
  }

  // TODO left shift binary operator>>(shift_by)

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

  static unsigned int getKmerSize()
  {
    return size;
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


protected:

  /**
   * @brief Sets all unused bits of the underlying k-mer data to 0.
   */
  inline void do_sanitize()
  {
    // TODO use templated helper struct for <0> template specialization
    data[nWords-1] &= getBitMask<word_type>(bitstream::invPadBits);
  }

  /**
   * @brief Sets all bits to zero.
   */
  inline void do_clear()
  {
    std::fill(data, data + nWords, static_cast<word_type>(0));
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
      this->do_left_shift(bitsPerChar);
      // copy `bitsperChar` least significant bits
      copyBitsFixed<word_type, bitsPerChar>(this->data[0], tmp_copy.data[0]);
      tmp_copy.do_right_shift(bitsPerChar);
    }

    // set ununsed bits to 0
    this->do_sanitize();
  }
};

} // namespace bliss

#endif // BLISS_COMMON_KMER_H
