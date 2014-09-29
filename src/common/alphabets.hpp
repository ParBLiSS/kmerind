/**
 * @file    alphabets.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares all common alphabets, including DNA, DNA5, RNA, RNA5
 *          AA (IUPAC), DNA_IUPAC, and CUSTOM
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_ALPHABETS_H
#define BLISS_COMMON_ALPHABETS_H

// own includes:
#include <common/base_types.hpp>

// TODO add the following alphabets:
// DNA
// RNA
// DNA5
// RNA5
// AA (IUPAC)
// DNA_IUPAC (15 characters)
// CUSTOM (no definition right now)

// need mapping from packed to unpacked.

// TODO add function calls for auto generation into documentation

struct BaseAlphabetChar
{
  // a castable char element for instaces of this struct
  CharType data_value;


  // TODO: assignment and construction are not directly inhertitable, thus those would need to be implemented in the separate alphabets
  //       (for general usage though, this would require virtual functions!????)

  // make this struct usable as a char
  operator CharType() const {return data_value;}


  /// copy assignment operators
  BaseAlphabetChar& operator=(const BaseAlphabetChar& c) {if (&c != this) data_value = c.data_value; return *this;}

  /// copy constructor
  BaseAlphabetChar(const BaseAlphabetChar& c) : data_value(c.data_value) {}

  // used by subclasses.
  BaseAlphabetChar& operator=(const CharType& c) {data_value = c; return *this;}
  BaseAlphabetChar(const CharType& c) : data_value(c) {}


  /// default constructor
  BaseAlphabetChar() {}


};


struct DNA : BaseAlphabetChar
{
  // TODO: check whether copying between char and DNA is actually optimized "away" (as it should)
  // This should make char and DNA useable interchangebly
  DNA& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
  DNA(const CharType& c) : BaseAlphabetChar(c) {}
  DNA() : BaseAlphabetChar() {}

  static constexpr AlphabetSizeType SIZE = 4;

  // lookup table for XYZ
  static constexpr uint8_t FROM_ASCII[256] =
  {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'A'     'C'             'G'
    0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,
//                 'T'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
//     'a'     'c'             'g'
    0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
//                 't'
    0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  };

  // reverse lookup table for XYZ
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T'  // = 3
  };

  // reverse lookup table for DNA5
  static constexpr char TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0  // = 3
  };


};

struct DNA5 : BaseAlphabetChar
{
    // TODO: check whether copying between char and DNA is actually optimized "away" (as it should)
    // This should make char and DNA useable interchangebly
    DNA5& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
    DNA5(const CharType& c) : BaseAlphabetChar(c) {}
    DNA5() : BaseAlphabetChar() {}

  static constexpr AlphabetSizeType SIZE = 5;
  static constexpr uint8_t FROM_ASCII[256] =
  {
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
//      'A'     'C'             'G'                         'N'
    4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
//                  'T'
    4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
//      'a'     'c'             'g'                         'n'
    4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
//                  't'
    4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
  };

  // reverse lookup table for DNA5
  static constexpr char TO_ASCII[SIZE] =
  {
    'A',  // = 0
    'C',  // = 1
    'G',  // = 2
    'T',  // = 3
    'N'  // = 4
  };

  // complement lookup table for DNA5
  static constexpr char TO_COMPLEMENT[SIZE] =
  {
    3,  // = 0
    2,  // = 1
    1,  // = 2
    0,  // = 3
    4   // = 4
  };
};

//constexpr uint8_t DNA::FROM_ASCII[256];
//constexpr char DNA::TO_ASCII[DNA::SIZE];
//constexpr uint8_t DNA5::FROM_ASCII[256];
//constexpr char DNA5::TO_ASCII[DNA5::SIZE];


#endif // BLISS_COMMON_ALPHABETS_H
