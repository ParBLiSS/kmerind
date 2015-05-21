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

// TODO add the following alphabets:
// D/RNA
// D/RNA5
// AA (IUPAC amino acid.)
// CUSTOM (no definition right now)

// TODO add function calls for auto generation into documentation
// TODO: check whether copying between char and DNA is actually optimized "away" (as it should)

namespace bliss
{
  namespace common
  {
    struct BaseAlphabetChar
    {
      // a castable char element for instaces of this struct
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
    
    ///  DNA alphabet: A T C G.  NOTE: DNA is configured so that negation is complement!
    struct DNA : BaseAlphabetChar
    {
      // This should make char and DNA useable interchangebly
      DNA& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
      DNA(const CharType& c) : BaseAlphabetChar(c) {}
      DNA() : BaseAlphabetChar() {}
    
      /// alphabet size
      static constexpr AlphabetSizeType SIZE = 4;
    
      /// ascii to alphabet lookup table
      static constexpr uint8_t FROM_ASCII[256] =
      {
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
      };
    
      /// alphabet to ascii lookup table
      static constexpr char TO_ASCII[SIZE] =
      {
        'A',  // = 0
        'C',  // = 1
        'G',  // = 2
        'T'  // = 3
      };
    
      /// complement lookup table
      static constexpr uint8_t TO_COMPLEMENT[SIZE] =
      {
        3,  // = 0    a->t
        2,  // = 1    c->g
        1,  // = 2    g->c
        0   // = 3    t->a
      };
    
      // this is same as reverse (by 2 bits) and invert for DNA, can do this in register. so don't need to use a lookup table.
      // 2 position in a nibble, 4 possible values each for 16 possible values.  going to bytes would exceed cacheline size
//      static constexpr uint8_t NIBBLE_REVERSE_COMPLEMENT[16] =
//      {
//       15,   // 0000, AA -> TT, 1111
//       11,   // 0001, AC -> GT, 1011
//        7,   // 0010, AG -> CT, 0111
//        3,   // 0011, AT -> AT, 0011
//       14,   // 0100, CA -> TG, 1110
//       10,   // 0101, CC -> GG, 1010
//        6,   // 0110, CG -> CG, 0110
//        2,   // 0111, CT -> AG, 0010
//       13,   // 1000, GA -> TC, 1101
//        9,   // 1001, GC -> GC, 1001
//        5,   // 1010, GG -> CC, 0101
//        1,   // 1011, GT -> AC, 0001
//       12,   // 1100, TA -> TA, 1100
//        8,   // 1101, TC -> GA, 1000
//        4,   // 1110, TG -> CA, 0100
//        0    // 1111, TT -> AA, 0000
//      };

      // do reverse in register so don't need to touch cache.
      // static constexpr uint8_t NIBBLE_REVERSE[16]


    };
    
    /// DNA5 Alphabet: A T C G N
    struct DNA5 : BaseAlphabetChar
    {
        // This should make char and DNA useable interchangebly
        DNA5& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        DNA5(const CharType& c) : BaseAlphabetChar(c) {}
        DNA5() : BaseAlphabetChar() {}
    
      /// alphabet size
      static constexpr AlphabetSizeType SIZE = 5;
    
      /// ascii to alphabet lookup table
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
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
      };
    
      /// alphabet to ascii lookup table
      static constexpr char TO_ASCII[SIZE] =
      {
        'A',  // = 0
        'C',  // = 1
        'G',  // = 2
        'T',  // = 3
        'N'  // = 4
      };
    
      /// complement lookup table
      static constexpr uint8_t TO_COMPLEMENT[SIZE] =
      {
        3,  // = 0
        2,  // = 1
        1,  // = 2
        0,  // = 3
        4   // = 4
      };

      // use default, serial reverse and rev complement

    };
    
    ///  DNA alphabet: A U C G.  NOTE: RNA is configured so that negation is complement!
    struct RNA : BaseAlphabetChar
    {
      // This should make char and RNA useable interchangebly
      RNA& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
      RNA(const CharType& c) : BaseAlphabetChar(c) {}
      RNA() : BaseAlphabetChar() {}
    
      /// alphabet size
      static constexpr AlphabetSizeType SIZE = 4;
    
      /// ascii to alphabet lookup table
      static constexpr uint8_t FROM_ASCII[256] =
      {
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
      };
    
      /// alphabet to ascii lookup table
      static constexpr char TO_ASCII[SIZE] =
      {
        'A',  // = 0
        'C',  // = 1
        'G',  // = 2
        'U'   // = 3
      };
    
      /// complement lookup table
      static constexpr uint8_t TO_COMPLEMENT[SIZE] =
      {
        3,  // = 0
        2,  // = 1
        1,  // = 2
        0   // = 3
      };
    

    };
    
    struct RNA5 : BaseAlphabetChar
    {
        // This should make char and RNA useable interchangebly
        RNA5& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        RNA5(const CharType& c) : BaseAlphabetChar(c) {}
        RNA5() : BaseAlphabetChar() {}
    
        /// alphabet size
      static constexpr AlphabetSizeType SIZE = 5;
    
      /// ascii to alphabet lookup table
      static constexpr uint8_t FROM_ASCII[256] =
      {
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
      };
    
      /// alphabet to ascii lookup table
      static constexpr char TO_ASCII[SIZE] =
      {
        'A',  // = 0
        'C',  // = 1
        'G',  // = 2
        'U',  // = 3
        'N'  // = 4
      };
    
      /// complement lookup table
      static constexpr uint8_t TO_COMPLEMENT[SIZE] =
      {
        3,  // = 0
        2,  // = 1
        1,  // = 2
        0,  // = 3
        4   // = 4
      };

      // use default, serial reverse and rev complement
    };
    
    
    
    /// IUPAC encoding, following strict interpretation of IUPAC-IUB SYMBOLS FOR NUCLEOTIDE (DNA OR RNA) NOMENCLATURE:
    ///                                                    Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030
    /// specifically, U is treated separately from T and X has same meaning as N.
    /// Note that because of U not being the same as T, the mapping between character and its complement is not bijective (no inverse).
    /// however, decoding preserves the original characters.   this is what BioPerl uses.
    struct DNA_IUPAC : BaseAlphabetChar
    {
        // TODO: check whether copying between char and DNA is actually optimized "away" (as it should)
        // This should make char and DNA useable interchangebly
        DNA_IUPAC& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
        DNA_IUPAC(const CharType& c) : BaseAlphabetChar(c) {}
        DNA_IUPAC() : BaseAlphabetChar() {}
    
        /// alphabet size
      static constexpr AlphabetSizeType SIZE = 16;
    
      /// ascii to alphabet lookup table
      static constexpr uint8_t FROM_ASCII[256] =
      {
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
      };
    
      /// alphabet to ascii lookup table
      /// note that ACGT each occupy 1 bit position.  U is specially represented as 0000.  all others are combinations of ACGT bits.
      static constexpr char TO_ASCII[SIZE] =
      {
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
      };
    
      /// complement lookup table
      static constexpr uint8_t TO_COMPLEMENT[SIZE] =
      {
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
      };
    };

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
      struct DNA16 : BaseAlphabetChar
      {
          // TODO: check whether copying between char and DNA is actually optimized "away" (as it should)
          // This should make char and DNA useable interchangebly
          DNA16& operator=(const CharType& c){ BaseAlphabetChar::operator=(c); return *this;}
          DNA16(const CharType& c) : BaseAlphabetChar(c) {}
          DNA16() : BaseAlphabetChar() {}

          /// alphabet size
        static constexpr AlphabetSizeType SIZE = 16;

        /// ascii to alphabet lookup table
        static constexpr uint8_t FROM_ASCII[256] =
        {
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
        };

        /// alphabet to ascii lookup table
        /// note that ACGT each occupy 1 bit position, indicating presence/absence.
        static constexpr char TO_ASCII[SIZE] =
        {
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
        };

        /// complement lookup table.  note that this is makes the complement a bit reversal.
        static constexpr uint8_t TO_COMPLEMENT[SIZE] =
        {
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
        };// good table of complements: http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html

    };
    
  } // namespace common
} // namespace bliss

#endif // BLISS_COMMON_ALPHABETS_H
