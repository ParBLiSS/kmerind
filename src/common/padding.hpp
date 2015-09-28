/**
 * @file    padding.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Implements common functions for managing padding.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_PADDING_H
#define BLISS_COMMON_PADDING_H

// C std lib includes:
#include <cstdlib>
#include <assert.h>

// C++ STL includes:
#include <iterator>
#include <limits>
#include <type_traits>

// own includes
#include "common/bit_ops.hpp"

namespace bliss
{
  namespace common
  {
  
    template <typename T, unsigned int BITS_PER_CHAR>
    struct PackingTraits
    {
      static_assert(std::is_integral<T>::value && !std::is_signed<T>::value, "Only supports unsigned integral types");
      /// The number of bits per data word (== number of bits when unsigned integral type)
      static constexpr unsigned int bits_per_word = std::numeric_limits<T>::digits;
      /// The number of bits per alphabet character (this is equal to the template
      /// parameter)
      static constexpr unsigned int bits_per_char = BITS_PER_CHAR;
      /// The number of data bits (bits per word minus padding bits) per storage word
      static constexpr unsigned int data_bits = roundDownToMultiple(bits_per_word, bits_per_char);
      /// The number of bits used as padding at the end of each data word.
      static constexpr unsigned int padding_bits = bits_per_word - data_bits;
      /// The number of characters inside of each storage/data word
      static constexpr unsigned int chars_per_word =  bits_per_word / bits_per_char;
    };
    
    /**
     * For stream of bits that is only padded at the very end of the stream
     * in case the number of total bits is not perfectly divisible by the
     * bits per word (of type T).
     */
    template <typename T, unsigned int NUM_TOTAL_BITS>
    struct UnpaddedStreamTraits
    {
        static_assert(std::is_integral<T>::value && !std::is_signed<T>::value, "Only supports unsigned integral types");
    public:
      /// The number of total bits in the stream
      static constexpr unsigned int nBits = NUM_TOTAL_BITS;
      /// The number of bytes needed to hold `nBits` bits
      static constexpr unsigned int nBytes = intCeil(nBits, 8);
      /// The number of bits for each word of type `T`
      static constexpr unsigned int bitsPerWord = std::numeric_limits<T>::digits;
      /// The number of bytes in each word of type `T`
      static constexpr unsigned int bytesPerWord = sizeof(T);
      /// The number of words of type `T` needed to hold `nBits` bits
      static constexpr unsigned int nWords = intCeil(nBytes, bytesPerWord);
      /// The number of padding bits at the end of the stream
      static constexpr unsigned int padBits = nWords*bitsPerWord - nBits;
      /// The number of data bits in the last word of the stream
      static constexpr unsigned int invPadBits = bitsPerWord - padBits;
      /// The number of bytes in the padding. This is the number of bytes that
      /// contain no data at all.
      static constexpr unsigned int padBytes = nWords*bytesPerWord - nBytes;
    };
    
    /*
    template <typename T, unsigned int NUM_TOTAL_BITS>
    constexpr unsigned int UnpaddedStreamTraits<T, NUM_TOTAL_BITS>::nBits;
    */
    
    /**
     * @brief Takes a sequence of integers (any size), removes padding in each word
     *        and concatenates the result into an output sequence of the same type.
     *
     * Note, that the output sequence is in general shorter than the input sequence,
     * because of padding bits removed.
     *
     * @tparam InputIterator  An input iterator (forward iterator concept).
     * @tparam OutputIterator An output iterator (forward iterator concept).
     *
     * @param begin[in|out]   The iterator to the begin of the input sequence.
     * @param out[in|out]     The iterator to the begin of the output sequence.
     * @param nBits[in]       The total number of bits that have to be read. This is
     *                        the termination condition of the function.
     * @param padBits[in]     The number of most significant bits used as padding in
     *                        each input value. This is the amount of bits removed
     *                        from each input word.
     * @param offset[in]      The bit offset for the current output word, to start
     *                        writing out bits. (default = 0)
     * @returns               The bit offset in the last word written (at position
     *                        `out`, all bits from this position up (towards the
     *                        most significant bits) are still unused.
     */
    template <typename InputIterator, typename OutputIterator>
    unsigned int removePaddingSameType(InputIterator& begin,
                                       OutputIterator& out,
                                       const unsigned int nBits,
                                       const unsigned int padBits,
                                       unsigned int offset = 0)
    {
      // get base types of iterators
      typedef typename std::iterator_traits<InputIterator>::value_type base_type;
      typedef typename std::iterator_traits<OutputIterator>::value_type out_type;
      // check that the size of the type is identical to the target type
      static_assert(sizeof(base_type) == sizeof(out_type),
                    "base and target type have to be of same byte size");
      // the offset has to be within the word size
      assert(offset < sizeof(base_type)*8);
    
      // the number of non padding bits per input word:
      const unsigned int bitsPerWord = sizeof(base_type)*8 - padBits;
      // the number of input words that have to be read to reach `nBits` bits
      unsigned int baseWordsToRead = intCeil(nBits, bitsPerWord);
      // the number of bits left to read
      unsigned int bitsToRead = nBits;
    
      // init the current offset
      unsigned int cur_offset = offset;
    
      // read all words from input sequence
      while (baseWordsToRead > 0)
      {
        // get the number of bits to read in this iteration:
        unsigned int readBits = std::min<unsigned int>(bitsPerWord, bitsToRead);
        // mask the next word with the number of bits to be read
        base_type nextWord = getLeastSignificantBitsMask<base_type>(readBits) & *begin;
    
        // write out the masked word
        if (cur_offset == 0)
          *out = nextWord;
        else
          // write out with offset
          *out |= nextWord << static_cast<base_type>(cur_offset);
    
        // check if the current word has overflown
        if (cur_offset + readBits >= sizeof(base_type)*8)
        {
          // update the offset for the next word
          cur_offset = (cur_offset + readBits) % (sizeof(base_type)*8);
          // increase output iterator
          ++out;
          // write out the cut off bits up to the next offset
          if (cur_offset > 0) *out = nextWord >> static_cast<base_type>(readBits - cur_offset);
        }
        else
        {
          // the current word still has empty bits left -> update offset
          cur_offset += readBits;
        }
    
        // update counters
        --baseWordsToRead;
        bitsToRead -= readBits;
        ++begin;
      }
    
      // return the next offset
      return cur_offset;
    }
    
    /**
     * @brief Removes padding of words in an input sequence and writes the
     *        concatenated sequence of bits into the output sequence.
     *
     * Note, that in this function, the base type of the input sequence and the
     * output sequence can be of different size.
     *
     * @tparam InputIterator  An input iterator type (forward iterator concept).
     * @tparam T              The base type of the output sequence.
     *
     * @param begin[in]   An iterator pointing to the first element in the input
     *                    sequence.
     * @param out[in]     A pointer pointing to the first element in the output
     *                    sequence.
     * @param nBits[in]   The number of bits the be read from the input sequence.
     * @param padBits[in] The number of padding bits (most significant bits) in
     *                    each input word.
     */
    template <typename InputIterator, typename T>
    void removePadding(InputIterator begin, T* out, const unsigned int nBits, const unsigned int padBits)
    {
      // TODO: change T* into an iterator by getting rid of all reinterpret_cast (see below)
      typedef typename std::iterator_traits<InputIterator>::value_type base_type;
      typedef T target_type;
    
      if (sizeof(target_type) >= sizeof(base_type))
      {
        // reinterpret the output pointer with the base type
        // TODO: this is cutting the input type into smaller pieces
        //       (could be implemented as a 1:2 or 1:4 (1:n) iterator
        base_type* target = reinterpret_cast<base_type*>(out);
        // simply copy all of them
        removePaddingSameType(begin, target, nBits, padBits);
      }
      else
      {
        const unsigned int baseBitsPerWord = sizeof(base_type)*8 - padBits;
        // the number of bits that can be read by the base type
        unsigned int baseBitsToRead = roundDownToMultiple(nBits, sizeof(base_type)*8);
        baseBitsToRead = roundDownToMultiple(baseBitsToRead, baseBitsPerWord);
    
    
        // TODO: do i need these?
        unsigned int remBits = nBits - baseBitsToRead;
    
        // TODO: this is combining the data into bigger pieces (n:1)
        base_type * target_base = reinterpret_cast<base_type*>(out);
    
        // copy all but the last word efficiently using the bigger base type
        unsigned int cur_offset = 0;
        if (baseBitsToRead > 0)
        {
          cur_offset = removePaddingSameType(begin, target_base, baseBitsToRead, padBits);
        }
    
        // the last word has to be handled differently
    
    
        base_type* base_to = target_base;
        // TODO: this is down casting (smaller type 1:n)
        target_type* to = reinterpret_cast<target_type*>(base_to);
        // skip already filled words
        unsigned int skip = cur_offset / (sizeof(target_type)*8);
        cur_offset %= sizeof(target_type)*8;
        to += skip;
    
        // iterate through all remaining words until all bits are read
        while (true)
        {
          // get next word
          base_type cur_base = *begin;
          // TODO: cast to pointer of the smaller target_type (1:n)
          target_type* from = reinterpret_cast<target_type*>(&cur_base);
    
          // get the number of bits the be read in the next step
          unsigned int nextNBits = std::min(remBits, baseBitsPerWord);
          // read the bits from the input pointer, with padding 0
          // TODO: this is a very special case (basically just copy a bunch of words
          //       over until nextNBits is met, which is the bits in the bigger word
          //       minus padding)
          cur_offset = removePaddingSameType(from, to, nextNBits, 0, cur_offset);
          if (remBits == nextNBits)
          {
            break;
          }
          else
          {
            remBits -= nextNBits;
            ++begin;
          }
        }
        //return remBits;
      }
    }
  
  
  } //namespace common
} // namespace bliss

#endif // BLISS_COMMON_PADDING_H
class A;
