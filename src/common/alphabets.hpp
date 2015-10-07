/**
 * @file    alphabets.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares all common alphabets, including DNA, DNA5, RNA, RNA5
 *          AA (IUPAC), DNA_IUPAC, and CUSTOM
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_ALPHABETS_H
#define BLISS_COMMON_ALPHABETS_H

// own includes:
#include "common/base_types.hpp"
#include "utils/constexpr_array.hpp"

// TODO add the following alphabets:
// D/RNA
// D/RNA5
// AA (IUPAC amino acid.)
// CUSTOM (no definition right now)

// TODO add function calls for auto generation into documentation
// TODO: check whether copying between char and DNA is actually optimized "away" (as it should)

#include <array>

namespace bliss
{
  namespace common
  {
    namespace alphabet
    {
    
      struct BaseAlphabetChar
      {
        // a castable char element for instances of this struct
        CharType data_value;

        // assignment and construction are not directly inhertitable, thus those would need to be implemented in the separate alphabets

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
    
      ///  ASCII alphabet.  note that TO_ASCII and FROM_ASCII do not do anything.  TO_COMPLEMENT more importantly does nothing.
      // dummy template parameter to take advantage of exemption for templated class - allowing definition of std::arrays to exist in header.  (declaration and initialization can be in header anyways)
      template <typename DUMMY = void>
      struct ASCII_T : BaseAlphabetChar
      {
        // This should make char and DNA useable interchangebly
        ASCII_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        ASCII_T(const CharType& c) : BaseAlphabetChar(c) {}
        ASCII_T() : BaseAlphabetChar() {}

        /// alphabet size
        static constexpr AlphabetSizeType SIZE = 256;

        /// ascii to alphabet lookup table
        static constexpr std::array<uint8_t, 256> FROM_ASCII = bliss::utils::make_array<uint8_t, 256>();

        /// alphabet to ascii lookup table
        static constexpr std::array<uint8_t, SIZE> TO_ASCII = bliss::utils::make_array<uint8_t, 256>();

        /// complement lookup table
        static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT = bliss::utils::make_array<uint8_t, 256>();

      };

      template <typename DUMMY>
      constexpr std::array<uint8_t, 256> ASCII_T<DUMMY>::FROM_ASCII;
      template <typename DUMMY>
      constexpr std::array<uint8_t, ASCII_T<DUMMY>::SIZE> ASCII_T<DUMMY>::TO_ASCII;
      template <typename DUMMY>
      constexpr std::array<uint8_t, ASCII_T<DUMMY>::SIZE> ASCII_T<DUMMY>::TO_COMPLEMENT;


      ///  DNA alphabet: A T C G.  NOTE: DNA is configured so that negation is complement!
      // dummy template parameter to take advantage of exemption for templated class - allowing definition of std::arrays to exist in header.  (declaration and initialization can be in header anyways)
      template <typename DUMMY = void>
      struct DNA_T : BaseAlphabetChar
      {
        // This should make char and DNA useable interchangebly
        DNA_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        DNA_T(const CharType& c) : BaseAlphabetChar(c) {}
        DNA_T() : BaseAlphabetChar() {}

        /// alphabet size
        static constexpr AlphabetSizeType SIZE = 4;

        /// ascii to alphabet lookup table
        static constexpr std::array<uint8_t, 256> FROM_ASCII =
        {{
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       //     'A'     'C'             'G'
           0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
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
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
         }};

        /// alphabet to ascii lookup table
        static constexpr std::array<char, SIZE> TO_ASCII =
        {{
          'A',  // = 0
          'C',  // = 1
          'G',  // = 2
          'T'  // = 3
        }};

        /// complement lookup table
        static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT =
        {{
          3,  // = 0    a->t
          2,  // = 1    c->g
          1,  // = 2    g->c
          0   // = 3    t->a
        }};
      };

      template <typename DUMMY>
      constexpr std::array<uint8_t, 256> DNA_T<DUMMY>::FROM_ASCII;
      template <typename DUMMY>
      constexpr std::array<char, DNA_T<DUMMY>::SIZE> DNA_T<DUMMY>::TO_ASCII;
      template <typename DUMMY>
      constexpr std::array<uint8_t, DNA_T<DUMMY>::SIZE> DNA_T<DUMMY>::TO_COMPLEMENT;



      /// DNA5 Alphabet: A T C G N
      template <typename DUMMY = void>
      struct DNA5_T : BaseAlphabetChar
      {
          // This should make char and DNA useable interchangebly
          DNA5_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
          DNA5_T(const CharType& c) : BaseAlphabetChar(c) {}
          DNA5_T() : BaseAlphabetChar() {}

        /// alphabet size
        static constexpr AlphabetSizeType SIZE = 5;

        /// ascii to alphabet lookup table
        static constexpr std::array<uint8_t, 256> FROM_ASCII =
        {{
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
          4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
        }};

        /// alphabet to ascii lookup table
        static constexpr std::array<char, SIZE> TO_ASCII =
        {{
          'A',  // = 0
          'C',  // = 1
          'G',  // = 2
          'T',  // = 3
          'N'  // = 4
        }};

        /// complement lookup table
        static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT =
        {{
          3,  // = 0
          2,  // = 1
          1,  // = 2
          0,  // = 3
          4   // = 4
        }};

      };
      template <typename DUMMY>
      constexpr std::array<uint8_t, 256> DNA5_T<DUMMY>::FROM_ASCII;
      template <typename DUMMY>
      constexpr std::array<char, DNA5_T<DUMMY>::SIZE> DNA5_T<DUMMY>::TO_ASCII;
      template <typename DUMMY>
      constexpr std::array<uint8_t, DNA5_T<DUMMY>::SIZE> DNA5_T<DUMMY>::TO_COMPLEMENT;



    ///  DNA alphabet: A U C G.  NOTE: RNA is configured so that negation is complement!
    template <typename DUMMY = void>
    struct RNA_T : BaseAlphabetChar
    {
      // This should make char and RNA useable interchangebly
      RNA_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
      RNA_T(const CharType& c) : BaseAlphabetChar(c) {}
      RNA_T() : BaseAlphabetChar() {}
    
      /// alphabet size
      static constexpr AlphabetSizeType SIZE = 4;
    
      /// ascii to alphabet lookup table
      static constexpr std::array<uint8_t, 256> FROM_ASCII =
      {{
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    //     'A'     'C'             'G'
        0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
    //                     'U'
        0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    //     'a'     'c'             'g'
        0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,
    //                     'u'
        0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
      }};
    
      /// alphabet to ascii lookup table
      static constexpr std::array<char, SIZE> TO_ASCII =
      {{
        'A',  // = 0
        'C',  // = 1
        'G',  // = 2
        'U'   // = 3
      }};
    
      /// complement lookup table
      static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT =
      {{
        3,  // = 0
        2,  // = 1
        1,  // = 2
        0   // = 3
      }};

    };
    
    template <typename DUMMY>
    constexpr std::array<uint8_t, 256> RNA_T<DUMMY>::FROM_ASCII;
    template <typename DUMMY>
    constexpr std::array<char, RNA_T<DUMMY>::SIZE> RNA_T<DUMMY>::TO_ASCII;
    template <typename DUMMY>
    constexpr std::array<uint8_t, RNA_T<DUMMY>::SIZE> RNA_T<DUMMY>::TO_COMPLEMENT;


    template <typename DUMMY = void>
    struct RNA5_T : BaseAlphabetChar
    {
        // This should make char and RNA useable interchangebly
        RNA5_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        RNA5_T(const CharType& c) : BaseAlphabetChar(c) {}
        RNA5_T() : BaseAlphabetChar() {}
    
        /// alphabet size
      static constexpr AlphabetSizeType SIZE = 5;
    
      /// ascii to alphabet lookup table
      static constexpr std::array<uint8_t, 256> FROM_ASCII =
      {{
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    //      'A'     'C'             'G'                         'N'
        4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
    //                     'U'
        4,  4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    //      'a'     'c'             'g'                         'n'
        4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
    //                     'u'
        4,  4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
      }};
    
      /// alphabet to ascii lookup table
      static constexpr std::array<char, SIZE> TO_ASCII =
      {{
        'A',  // = 0
        'C',  // = 1
        'G',  // = 2
        'U',  // = 3
        'N'  // = 4
      }};
    
      /// complement lookup table
      static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT =
      {{
        3,  // = 0
        2,  // = 1
        1,  // = 2
        0,  // = 3
        4   // = 4
      }};

      // use default, serial reverse and rev complement
      static constexpr uint8_t from_ascii(uint8_t ascii) {
        return FROM_ASCII[ascii];
      }
    };
    
    template <typename DUMMY>
    constexpr std::array<uint8_t, 256> RNA5_T<DUMMY>::FROM_ASCII;
    template <typename DUMMY>
    constexpr std::array<char, RNA5_T<DUMMY>::SIZE> RNA5_T<DUMMY>::TO_ASCII;
    template <typename DUMMY>
    constexpr std::array<uint8_t, RNA5_T<DUMMY>::SIZE> RNA5_T<DUMMY>::TO_COMPLEMENT;

    
    /// IUPAC encoding, following strict interpretation of IUPAC-IUB SYMBOLS FOR NUCLEOTIDE (DNA OR RNA) NOMENCLATURE:
    ///                                                    Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030
    /// specifically, U is treated separately from T and X has same meaning as N.
    /// Note that because of U not being the same as T, the mapping between character and its complement is not bijective (no inverse).
    /// however, decoding preserves the original characters.   this is what BioPerl uses.
    template <typename DUMMY = void>
    struct DNA_IUPAC_T : BaseAlphabetChar
    {
        // TODO: check whether copying between char and DNA is actually optimized "away" (as it should)
        // This should make char and DNA useable interchangebly
        DNA_IUPAC_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        DNA_IUPAC_T(const CharType& c) : BaseAlphabetChar(c) {}
        DNA_IUPAC_T() : BaseAlphabetChar() {}
    
        /// alphabet size
      static constexpr AlphabetSizeType SIZE = 16;
    
      /// ascii to alphabet lookup table
      static constexpr std::array<uint8_t, 256> FROM_ASCII =
      {{
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        //   'A'  'B'  'C'  'D'            'G'  'H'            'K'       'M'  'N'
        0xF, 0x1, 0xE, 0x2, 0xD, 0xF, 0xF, 0x4, 0xB, 0xF, 0xF, 0xC, 0xF, 0x3, 0xF, 0xF,
        //        'R'  'S'  'T'  'U'  'V'  'W'  'X'  'Y'
        0xF, 0xF, 0x5, 0x6, 0x8, 0x0, 0x7, 0x9, 0xF, 0xA, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        //   'a'  'b'  'c'  'd'            'g'  'h'            'k'       'm'  'n'
        0xF, 0x1, 0xE, 0x2, 0xD, 0xF, 0xF, 0x4, 0xB, 0xF, 0xF, 0xC, 0xF, 0x3, 0xF, 0xF,
        //        'r'  's'  't'  'u'  'v'  'w'  'x'  'y'
        0xF, 0xF, 0x5, 0x6, 0x8, 0x0, 0x7, 0x9, 0xF, 0xA, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
        0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF
      }};
    
      /// alphabet to ascii lookup table
      /// note that ACGT each occupy 1 bit position.  U is specially represented as 0000.  all others are combinations of ACGT bits.
      static constexpr std::array<char, SIZE> TO_ASCII =
      {{
        'U',  // = 0      0000
        'A',  // = 1      0001
        'C',  // = 2      0010
        'M',  // = 3      0011
        'G',  // = 4      0100
        'R',  // = 5      0101
        'S',  // = 6      0110
        'V',  // = 7      0111
        'T',  // = 8      1000
        'W',  // = 9      1001
        'Y',  // = 10     1010
        'H',  // = 11     1011
        'K',  // = 12     1100
        'D',  // = 13     1101
        'B',  // = 14     1110
        'N'   // = 15     1111
      }};
    
      /// complement lookup table
      static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT =
      {{
        1,  // = 0  0000  U      U->A  A            0001
        8,  // = 1  0001  A      A->T  T    (no TU) 1000
        4,  // = 2  0010  C      C->G  G            0100
        12, // = 3  0011  AC     M->K  TUG          1100
        2,  // = 4  0100  G      G->C  C            0010
        10, // = 5  0101  AG     R->Y  CTU          1010
        6,  // = 6  0110  CG     S->S  CG           0110
        14, // = 7  0111  ~T~U   V->B  ~A           1110
        1,  // = 8  1000  T      T->A  A            0001
        9,  // = 9  1001  ATU    W->W  ATU          1001
        5,  // = 10 1010  CTU    Y->R  AG           0101
        13, // = 11 1011  ~G     H->D  ~C           1101
        3,  // = 12 1100  TUG    K->M  AC           0011
        11, // = 13 1101  ~C     D->H  ~G           1011
        7,  // = 14 1110  ~A     B->V  ~T~U         0111
        15  // = 15 1111  N      N->N  N            1111
      }};

    };


    template <typename DUMMY>
    constexpr std::array<uint8_t, 256> DNA_IUPAC_T<DUMMY>::FROM_ASCII;
    template <typename DUMMY>
    constexpr std::array<char, DNA_IUPAC_T<DUMMY>::SIZE> DNA_IUPAC_T<DUMMY>::TO_ASCII;
    template <typename DUMMY>
    constexpr std::array<uint8_t, DNA_IUPAC_T<DUMMY>::SIZE> DNA_IUPAC_T<DUMMY>::TO_COMPLEMENT;



      /// Looser IUPAC encoding
      ///
      /// specifically, U == T, X == N, and a gap character ('.' or '-') is allowed, which makes it more friendly for alignment.
      ///
      /// note that the mapping between a character and its complement is bijective and invertible.
      /// however, decoding will NOT return the original string when U and T are involved.
      ///
      /// The encoding has been choosen so that bit reversal creates complement.
      ///
      /// This is what is used by bioinformatics.org (T==U), dnabaser.com ('-' == gap),
      /// wikipedia (reporting '-' == N), and finally and definitively,
      /// http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html#301 (15 chars only, no U) (NC-IUB recommendation 1984)
      template <typename DUMMY = void>
      struct DNA16_T : BaseAlphabetChar
      {
          // TODO: check whether copying between char and DNA is actually optimized "away" (as it should)
          // This should make char and DNA useable interchangebly
          DNA16_T& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
          DNA16_T(const CharType& c) : BaseAlphabetChar(c) {}
          DNA16_T() : BaseAlphabetChar() {}

          /// alphabet size
        static constexpr AlphabetSizeType SIZE = 16;

        /// ascii to alphabet lookup table
        static constexpr std::array<uint8_t, 256> FROM_ASCII =
        {{
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
       // ' '                                                              '-'  '.'
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0x0, 0x0, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          //   'A'  'B'  'C'  'D'            'G'  'H'            'K'       'M'  'N'
          0xF, 0x1, 0xE, 0x2, 0xD, 0xF, 0xF, 0x4, 0xB, 0xF, 0xF, 0xC, 0xF, 0x3, 0xF, 0xF,
          //        'R'  'S'  'T'  'U'  'V'  'W'  'X'  'Y'
          0xF, 0xF, 0x5, 0x6, 0x8, 0x8, 0x7, 0x9, 0xF, 0xA, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          //   'a'  'b'  'c'  'd'            'g'  'h'            'k'       'm'  'n'
          0xF, 0x1, 0xE, 0x2, 0xD, 0xF, 0xF, 0x4, 0xB, 0xF, 0xF, 0xC, 0xF, 0x3, 0xF, 0xF,
          //        'r'  's'  't'  'u'  'v'  'w'  'x'  'y'
          0xF, 0xF, 0x5, 0x6, 0x8, 0x8, 0x7, 0x9, 0xF, 0xA, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
          0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF
        }};

        /// alphabet to ascii lookup table
        /// note that ACGT each occupy 1 bit position, indicating presence/absence.
        static constexpr std::array<char, SIZE> TO_ASCII =
        {{
          '.',  // = 0   0x0   0000  // gap character allowed.
          'A',  // = 1   0x1   0001
          'C',  // = 2   0x2   0010
          'M',  // = 3   0x3   0011
          'G',  // = 4   0x4   0100
          'R',  // = 5   0x5   0101
          'S',  // = 6   0x6   0110
          'V',  // = 7   0x7   0111
          'T',  // = 8   0x8   1000   // choose to use T instead of U
          'W',  // = 9   0x9   1001
          'Y',  // = 10  0xA   1010
          'H',  // = 11  0xB   1011
          'K',  // = 12  0xC   1100
          'D',  // = 13  0xD   1101
          'B',  // = 14  0xE   1110
          'N'   // = 15  0xF   1111   // choose to use N instead of X.
        }};

        /// complement lookup table.  note that this is makes the complement a bit reversal.
        static constexpr std::array<uint8_t, SIZE> TO_COMPLEMENT =
        {{
          0,  // = 0  0000  .      .->.  gap          0000
          8,  // = 1  0001  A      A->T  T    (no TU) 1000
          4,  // = 2  0010  C      C->G  G            0100
          12, // = 3  0011  AC     M->K  TUG          1100
          2,  // = 4  0100  G      G->C  C            0010
          10, // = 5  0101  AG     R->Y  CTU          1010
          6,  // = 6  0110  CG     S->S  CG           0110
          14, // = 7  0111  ~T~U   V->B  ~A           1110
          1,  // = 8  1000  T      T->A  A            0001
          9,  // = 9  1001  ATU    W->W  ATU          1001
          5,  // = 10 1010  CTU    Y->R  AG           0101
          13, // = 11 1011  ~G     H->D  ~C           1101
          3,  // = 12 1100  TUG    K->M  AC           0011
          11, // = 13 1101  ~C     D->H  ~G           1011
          7,  // = 14 1110  ~A     B->V  ~T~U         0111
          15  // = 15 1111  N      N->N  N            1111
        }};// good table of complements: http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html


    };
      template <typename DUMMY>
      constexpr std::array<uint8_t, 256> DNA16_T<DUMMY>::FROM_ASCII;
      template <typename DUMMY>
      constexpr std::array<char, DNA16_T<DUMMY>::SIZE> DNA16_T<DUMMY>::TO_ASCII;
      template <typename DUMMY>
      constexpr std::array<uint8_t, DNA16_T<DUMMY>::SIZE> DNA16_T<DUMMY>::TO_COMPLEMENT;


    } // namespace alphabet

      using ASCII = ::bliss::common::alphabet::ASCII_T<>;
      using DNA = ::bliss::common::alphabet::DNA_T<>;
      using DNA5 = ::bliss::common::alphabet::DNA5_T<>;
      using RNA = ::bliss::common::alphabet::RNA_T<>;
      using RNA5 = ::bliss::common::alphabet::RNA5_T<>;
      using DNA16 = ::bliss::common::alphabet::DNA16_T<>;
      using DNA_IUPAC = ::bliss::common::alphabet::DNA_IUPAC_T<>;


  } // namespace common
} // namespace bliss

#endif // BLISS_COMMON_ALPHABETS_H
