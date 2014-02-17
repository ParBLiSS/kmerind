/**
 * @file    PackingIterator.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_PACKINGITERATOR_H
#define BLISS_COMMON_PACKINGITERATOR_H

// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <iterator>

// own includes
#include <common/base_types.hpp>
#include <common/bit_ops.hpp>


namespace bliss
{


// TODO add template specialization for padding (TRUE or FALSE)
template <typename BaseIterator, BitSizeType bits_per_char, typename PackedStorageType=WordType>
class PackingIterator
  : public std::iterator<std::forward_iterator_tag,//typename std::iterator_traits<BaseIterator>::iterator_category,
                  PackedStorageType,
                  typename std::iterator_traits<BaseIterator>::difference_type>
{
public:
  typedef PackedStorageType value_type;
  typedef BaseIterator base_iterator;
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;
  typedef PackingIterator<BaseIterator, bits_per_char, PackedStorageType> this_type;

  // bits per character in the packed representation
  static constexpr BitSizeType bitsPerChar = bits_per_char;
  // chars that fit in each word
  static constexpr BitSizeType charsPerWord = (sizeof(PackedStorageType) * 8) / bitsPerChar;
  // how much padding is added to each word
  static constexpr BitSizeType paddingPerWord = sizeof(PackedStorageType) * 8 - bitsPerChar*charsPerWord;



  // default contructor:
  PackingIterator() {}

  PackingIterator(BaseIterator baseBegin)
    : baseCur(baseBegin), baseNext(baseBegin), baseEnd(baseBegin), buffer(0) {}

  PackingIterator(BaseIterator baseBegin, BaseIterator baseEnd)
    : baseCur(baseBegin), baseNext(baseBegin), baseEnd(baseEnd),
      buffer(0) {}


  const value_type& operator*()
  {
    // return packed word
    if (baseCur == baseNext)
    {
      // advance the next iterator and buffer result
      BitSizeType chars = charsPerWord;
      BitSizeType numChars = static_cast<BitSizeType>(std::min(static_cast<size_t>(chars), static_cast<size_t>(std::distance(baseNext, baseEnd))));
      std::array<base_value_type, charsPerWord> word_cache;
      for (BitSizeType i = 0; i < numChars; ++i)
      {
        word_cache[i] = *baseNext;
        ++baseNext;
      }
      // read them in reverse
      buffer = 0;
      for (BitSizeType i = numChars; i > 0; --i)
      //for (BitSizeType i = 1; i <= numChars; ++i)
      {
        buffer <<= bitsPerChar;
        buffer |= getBitMask<PackedStorageType>(bitsPerChar) & word_cache[i-1];
      }
    }
    return buffer;
  }

  this_type& operator++() // prefix ++
  {
    // in case the iterator has not been read yet, advance the next iterator
    if (baseCur == baseNext)
      // FIXME: replace advance(min()) with a own advance that increases while != baseEnd
      //        otherwise this might be in O(nk) instead of O(n)
      std::advance(baseNext, std::min(std::distance(baseNext, baseEnd), static_cast<typename std::iterator_traits<BaseIterator>::difference_type>(charsPerWord)));
    baseCur = baseNext;
    return (*this);
  }

  // TODO: implement operator--(), operator-(diff), operator+(diff)

  bool operator==(const this_type& rhs)
  {
    return this->baseCur == rhs.baseCur;
  }

  bool operator!=(const this_type& rhs)
  {
    return !this->operator==(rhs);
  }

  this_type& operator++(int) // postfix ++
  {
    // create a copy prior to increasing
    this_type tmp(*this);
    this->operator++();
    return tmp;
  }

private:

  // the base iterator
  BaseIterator baseCur;

  BaseIterator baseNext;

  BaseIterator baseEnd;

  value_type buffer;
};



// TODO add template specialization for padding (TRUE or FALSE)
template <typename BaseIterator, BitSizeType bits_per_char, typename UnpackedStorageType=CharType>
class UnpackingIterator
  : public std::iterator<std::forward_iterator_tag,//typename std::iterator_traits<BaseIterator>::iterator_category,
                  UnpackedStorageType,
                  typename std::iterator_traits<BaseIterator>::difference_type>
{
public:
  UnpackingIterator(BaseIterator baseBegin)
    : baseCur(baseBegin), curCharOffset(0) {}
  UnpackingIterator(BaseIterator baseBegin, unsigned int offset)
  {
    unsigned int nWords = offset / charsPerWord;
    curCharOffset = offset % charsPerWord;
    baseCur = baseBegin;
    std::advance(baseCur, nWords);
  }

  typedef UnpackedStorageType value_type;
  typedef BaseIterator base_iterator;
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;
  typedef UnpackingIterator<BaseIterator, bits_per_char, UnpackedStorageType> this_type;


  value_type operator*()
  {
    base_value_type x = getBitMask<base_value_type>(bitsPerChar) & (*baseCur >> (curCharOffset*bitsPerChar));
    return static_cast<value_type>(x);
  }

  this_type& operator++() // prefix ++
  {
    // increment
    ++curCharOffset;
    if (curCharOffset == charsPerWord)
    {
      ++baseCur;
      curCharOffset = 0;
    }
    return (*this);
  }

  // TODO: implement operator--(), operator-(diff), operator+(diff)

  bool operator==(const this_type& rhs)
  {
    return (this->baseCur == rhs.baseCur) && (this->curCharOffset == rhs.curCharOffset);
  }

  bool operator!=(const this_type& rhs)
  {
    return !this->operator==(rhs);
  }

  this_type& operator++(int) // postfix ++
  {
    // create a copy prior to increasing
    this_type tmp(*this);
    this->operator++();
    return tmp;
  }

private:
  // TODO maybe unify these into one shared super class
  // bits per character in the packed representation
  static constexpr BitSizeType bitsPerChar = bits_per_char;
  // chars that fit in each word (with padding)
  static constexpr BitSizeType charsPerWord = (sizeof(base_value_type) * 8) / bitsPerChar;
   // how much padding is added to each word
  static constexpr BitSizeType paddingPerWord = sizeof(base_value_type) * 8 - bitsPerChar*charsPerWord;

  // the base iterator
  BaseIterator baseCur;

  // the current offset in chars within one storage word
  BitSizeType curCharOffset = 0;
};

} // namespace bliss

#endif // BLISS_COMMON_PACKINGITERATOR_H