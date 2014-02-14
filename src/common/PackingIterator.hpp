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
  PackingIterator(BaseIterator baseBegin, BaseIterator baseEnd)
    : baseCur(baseBegin), baseNext(baseBegin), baseEnd(baseEnd),
      buffer(0) {}

  typedef PackedStorageType value_type;

  const value_type& operator*()
  {
    // return packed word
    if (baseCur == baseNext)
    {
      buffer = 0;
      // advance the next iterator and buffer result
      for (BitSizeType i = 0; i < charsPerWord && baseNext != baseEnd; ++i)
      {
        buffer <<= bitsPerChar;
        buffer |= getBitMask<PackedStorageType>(bitsPerChar) & *baseNext;
        ++baseNext;
      }
    }
    return buffer;
  }

  PackingIterator<BaseIterator, bits_per_char, PackedStorageType>& operator++() // prefix ++
  {
    // in case the iterator has not been read yet, advance the next iterator
    if (baseCur == baseNext)
      std::advance(baseNext, std::min(std::distance(baseNext, baseEnd), static_cast<typename std::iterator_traits<BaseIterator>::difference_type>(charsPerWord)));
    baseCur = baseNext;
    return (*this);
  }

  // TODO: implement operator--(), operator-(diff), operator+(diff)

  bool operator==(const PackingIterator<BaseIterator, bits_per_char, PackedStorageType>& rhs)
  {
    return this->baseCur == rhs.baseCur;
  }

  bool operator!=(const PackingIterator<BaseIterator, bits_per_char, PackedStorageType>& rhs)
  {
    return !this->operator==(rhs);
  }

  PackingIterator<BaseIterator, bits_per_char, PackedStorageType>& operator++(int) // postfix ++
  {
    // create a copy prior to increasing
    PackingIterator<BaseIterator, bits_per_char, PackedStorageType> tmp(*this);
    this->operator++();
    return tmp;
  }

private:
  // bits per character in the packed representation
  static constexpr BitSizeType bitsPerChar = bits_per_char;
  // chars that fit in each word
  static constexpr BitSizeType charsPerWord = (sizeof(PackedStorageType) * 8) / bitsPerChar;
  // how much padding is added to each word
  static constexpr BitSizeType paddingPerWord = sizeof(PackedStorageType) * 8 - bitsPerChar*charsPerWord;

  // the base iterator
  BaseIterator baseCur;

  BaseIterator baseNext;

  BaseIterator baseEnd;

  value_type buffer;
};

} // namespace bliss

#endif // BLISS_COMMON_PACKINGITERATOR_H