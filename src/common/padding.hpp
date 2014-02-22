/**
 * @file    padding.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief
 *
 * Copyright (c) TODO
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

// own includes
#include <common/bit_ops.hpp>

namespace bliss
{

template <typename InputIterator, typename OutputIterator>
unsigned int removePaddingCopyHelper(InputIterator& begin, OutputIterator& out, const unsigned int nBits, const unsigned int padBits, unsigned int offset = 0)
{
  // get base types of iterators
  typedef typename std::iterator_traits<InputIterator>::value_type base_type;
  typedef typename std::iterator_traits<OutputIterator>::value_type out_type;
  // check that the size of the type is identical to the target type
  static_assert(sizeof(base_type) == sizeof(out_type), "base and target type have to be of same byte size");
  // the offset has to be within the word size
  assert(offset < sizeof(base_type)*8);

  // the number of non padding bits per input word:
  const unsigned int bitsPerWord = sizeof(base_type)*8 - padBits;
  // the number of input words that have to be read to reach `nBits` bits
  unsigned int baseWordsToRead = intCeil<unsigned int>(nBits, bitsPerWord);
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
    base_type nextWord = getBitMask<base_type>(readBits) & *begin;

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
      *out = nextWord >> static_cast<base_type>(readBits - cur_offset);
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


template <typename InputIterator, typename T>
void copyRemPadding(InputIterator begin, T* out, const unsigned int nBits, const unsigned int padBits)
{
  typedef typename std::iterator_traits<InputIterator>::value_type base_type;
  typedef T target_type;

  if (sizeof(target_type) >= sizeof(base_type))
  {
    // reinterpret the output pointer with the base type
    // TODO: this is cutting the input type into smaller pieces
    //       (could be implemented as a 1:2 or 1:4 (1:n) iterator
    base_type* target = reinterpret_cast<base_type*>(out);
    // simply copy all of them
    removePaddingCopyHelper(begin, target, nBits, padBits);
  }
  else
  {
    const unsigned int baseBitsPerWord = sizeof(base_type)*8 - padBits;
    // the number of bits that can be read by the base type
    unsigned int baseBitsToRead = roundDownToMultiple<unsigned int>(nBits, sizeof(base_type)*8);
    baseBitsToRead = roundDownToMultiple<unsigned int>(baseBitsToRead, baseBitsPerWord);


    // TODO: do i need these?
    unsigned int remBits = nBits - baseBitsToRead;

    // TODO: this is combining the data into bigger pieces (n:1)
    base_type * target_base = reinterpret_cast<base_type*>(out);

    // copy all but the last word
    unsigned int cur_offset = 0;
    if (baseBitsToRead > 0)
    {
      cur_offset = removePaddingCopyHelper(begin, target_base, baseBitsToRead, padBits);
    }

    // the last word has to be handled differently


    base_type* base_to = target_base;
    // TODO: this is down casting (smaller type 1:n)
    target_type* to = reinterpret_cast<target_type*>(base_to);
    unsigned int skip = cur_offset / (sizeof(target_type)*8);
    cur_offset %= sizeof(target_type)*8;
    to += skip;

    while (remBits > 0)
    {
      base_type cur_base = *begin;
      // TODO: cast to the smaller target_type (1:n)
      target_type* from = reinterpret_cast<target_type*>(&cur_base);


      unsigned int nextNBits = std::min(remBits, baseBitsPerWord);
      cur_offset = removePaddingCopyHelper(from, to, nextNBits, 0, cur_offset);
      remBits -= nextNBits;
      ++begin;
    }
  }
}


} // namespace bliss

#endif // BLISS_COMMON_PADDING_H
class A;