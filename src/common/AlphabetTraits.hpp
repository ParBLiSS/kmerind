/**
 * @file    AlphabetTraits.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares the class `AlphabetTraits<T>`, which implements static
 *          functions returning the traits of alphabet `T`.
 *          These traits include: the size of the alphabet, the bitsize
 *          per character of the alphabet, the PackedString type for the
 *          alphabet, and translation functions for the alphabet (translating
 *          between ASCII and internal representation).
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_ALPHABETTRAITS_H
#define BLISS_COMMON_ALPHABETTRAITS_H

// C std lib includes
#include <cstdlib>

// own includes:
#include <common/bit_ops.hpp>
#include <common/alphabets.hpp>
#include <common/PackedString.hpp>

template <typename T>
struct AlphabetTraits {
  // TODO: implement a static_assert to check for all of the needed members of T
  // http://stackoverflow.com/questions/257288/is-it-possible-to-write-a-c-template-to-check-for-a-functions-existence
  // (i.e. to ensure everything is properly implemented)
  static_assert(sizeof(T::SIZE), "the alphabet T must implement the constexpr `SIZE`"); 
  static_assert(sizeof(T::FROM_ASCII)/sizeof(T::FROM_ASCII[0]) == 256, "the alphabet T must provide a translation table `FROM_ASCII` for all 256 ASCII values.");
  static_assert(sizeof(T::TO_ASCII)/sizeof(T::TO_ASCII[0]) == T::SIZE, "the alphabet T must provide a translation table `TO_ASCII` of the same length as `SIZE`.");
  static_assert(sizeof(T) == 1, "The alphabet base class cannot be of size != 1");

private:
  static constexpr AlphabetSizeType SIZE = T::SIZE;
  static constexpr BitSizeType BITS_PER_CHAR = bits_per_char_needed(T::SIZE);

public:
  typedef PackedStringImpl<BITS_PER_CHAR> PackedStringType;

  static constexpr unsigned int getSize()
  {
    return SIZE;
  }

  static constexpr unsigned int getBitsPerChar()
  {
    return BITS_PER_CHAR;
  }

  template<typename InputIterator, typename OutputIterator>
  static OutputIterator translateFromAscii(InputIterator first, InputIterator last, OutputIterator result)
  {
    // TODO implement static_assert's to check for iterator properties
    while (first != last) {
      *(result++) = T::FROM_ASCII[static_cast<size_t>(*(first++))];
    }
    return result;
  }

  template<typename InputIterator, typename OutputIterator>
  static OutputIterator translateToAscii(InputIterator first, InputIterator last, OutputIterator result)
  {
    // TODO implement static_assert's to check for iterator properties
    while (first != last) {
      // TODO proper ASSERT handleing
      //assert(*first < SIZE, "the alphabet value is out of range for back-translation to ASCII");
      *(result++) = T::TO_ASCII[static_cast<size_t>(*(first++))];
    }
    return result;
  }
};


#endif // BLISS_COMMON_ALPHABETTRAITS_H
