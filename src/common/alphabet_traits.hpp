/**
 * @file    alphabet_traits.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief   Declares the class `AlphabetTraits<T>`, which implements static
 *          functions returning the traits of alphabet `T`.
 *          These traits include: the size of the alphabet, the bitsize
 *          per character of the alphabet, the PackedString type for the
 *          alphabet, and translation functions for the alphabet (translating
 *          between ASCII and internal representation).
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_ALPHABETTRAITS_H
#define BLISS_COMMON_ALPHABETTRAITS_H

// C std lib includes
#include <cstdlib>

// own includes:
#include "common/bit_ops.hpp"
#include "common/alphabets.hpp"
#include "common/packed_string.hpp"


namespace bliss
{
  namespace common
  {
    /// Functor to convert from ascii to values in the alphabet
    template<typename Alphabet, typename I = unsigned char>
    struct ASCII2 {
        uint8_t operator()(I ascii) const {
          return Alphabet::FROM_ASCII[static_cast<size_t>(ascii)];
        }
    };
  
    /// Functor to convert from values in the alphabet to is complement.
    template<typename Alphabet>
    struct ToComplement {
        uint8_t operator()(uint8_t in) const {
          return Alphabet::TO_COMPlEMENT[static_cast<size_t>(in)];
        }
    };
  
  /**
   * @brief Templated class of static functions returning properties of the
   *        alphabet.
   *
   * This class is templated by the alphabet and will provide the properties
   * of the alphabet via it's static member functions. The properties include
   * the alphabet's size, the bitsize per character. Furthermore, this class
   * supplies methods for translating from and to ASCII representation.
   * An example use:
   * @code
   *    std::basic_string<DNA5> read = "ACTGCATGGAAC"; // TODO: this might not work
   *    // getting some basic information:
   *    int alphabet_size = AlphabetTraits<DNA5>::getSize(); // will return 5
   *    int bits = AlphabetTraits<DNA5>::getBitsPerChar(); // will return 3 (bits)
   *    // translating from ascii into internal representation (needed before
   *    // packing of the sequence into a PackedString)
   *    AlphabetTraits<DNA5>::translateFromAscii(read.begin(), read.end(), read.begin());
   * @endcode
   *
   */
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
    /// Holds the size of the alphabet
    static constexpr AlphabetSizeType SIZE = T::SIZE;
    /// Holds the bit size of the alphabet, this is the number of bits that is
    /// needed to represent each symbol of the alphabet
    static constexpr BitSizeType BITS_PER_CHAR = ceilLog2(T::SIZE);
  
    // private contructor: don't allow initialization
    AlphabetTraits() {};
  
  public:
    /// The type for the PackedString compatible with this alphabet
    typedef PackedStringImpl<BITS_PER_CHAR> PackedStringType;
  
  
    /**
     * @brief  Returns the size of the alphabet.
     *
     * The size of the alphabet is defined to be the number of unique symbols of
     * the alphabet.
     *
     * @return The size of the alphabet.
     */
    static constexpr unsigned int getSize()
    {
      return SIZE;
    }
  
    /**
     * @brief Returns the bits per character of the alphabet.
     *
     * The number of bits needed to represent each character of the alphabet is
     * given by the ceiling of the log2 of the size of the alphabet. This is
     * exactly what this function returns.
     *
     * @return The bits per character of the alphabet.
     */
    static constexpr unsigned int getBitsPerChar()
    {
      return BITS_PER_CHAR;
    }
  
  
    /**
     * @brief Translates the given input sequence from ASCII format to internal
     *        format used by the alphabet.
     *
     * The ASCII characters are mapped to values in the range `{0,...,size -1}`,
     * where `size` is the size of the alphabet. This function translates the
     * input in-place.
     *
     * @param first   InputIterator to the start of the sequence to be translated.
     * @param last    InputIterator to the end of the sequence to be translated
     *                (that is one element past the last element in the sequence).
     *
     * @tparam InputOutputIterator  Iterator type for the `begin()` and `end()`
     *                              iterators.
     */
    template<typename InputIterator, typename OutputIterator>
    static void translateFromAscii(InputIterator first, InputIterator last)
    {
      translateFromAscii(first, last, first);
    }
  
    /**
     * @brief Translates the given input sequence from ASCII format to internal
     *        format used by the alphabet.
     *
     * The ASCII characters are mapped to values in the range `{0,...,size -1}`,
     * where `size` is the size of the alphabet. This function also works
     * in-place, i.e. `first` can be equal to `result`. In the case of overlapping
     * sequences with `first` != `result`, the behavious is not defined.
     *
     * @param first   InputIterator to the start of the sequence to be translated.
     * @param last    InputIterator to the end of the sequence to be translated
     *                (that is one element past the last element in the sequence).
     * @param result  OuptutIterator to the start of the result sequence.
     *
     * @tparam InputIterator  InputIterator type for the `begin()` and `end()`
     *                        iterators.
     * @tparam OutputIterator Iterator type for sequential output operations.
     *
     * @return        An Iterator to the end of the result sequence, this will
     *                point to one element past the last written element.
     */
    template<typename InputIterator, typename OutputIterator>
    static OutputIterator translateFromAscii(InputIterator first, InputIterator last, OutputIterator result)
    {
      // TODO implement static_assert's to check for iterator properties
      while (first != last) {
        *(result++) = T::FROM_ASCII[static_cast<size_t>(*(first++))];
      }
      return result;
    }
  
    /**
     * @brief Translates the given input sequence from the internal format to
     *        ASCII format.
     *
     * Values in the range `{0,...,size -1}` (where `size` is the size of the
     * alphabet) are mapped to ASCII characters, as is defined by the alphabet.
     * E.g. for DNA, this will map {0 -> 'A', 1 -> 'C', ...}.
     * This function translates the given sequence in-place.
     *
     * @param first   InputIterator to the start of the sequence to be translated.
     * @param last    InputIterator to the end of the sequence to be translated
     *                (that is one element past the last element in the sequence).
     *
     * @tparam InputIterator  Iterator type for the `begin()` and `end()`
     *                        iterators.
     */
    template<typename InputIterator>
    static void translateToAscii(InputIterator first, InputIterator last)
    {
      translateToAscii(first, last, first);
    }
  
    /**
     * @brief Translates the given input sequence from the internal format to
     *        ASCII format.
     *
     * Values in the range `{0,...,size -1}` (where `size` is the size of the
     * alphabet) are mapped to ASCII characters, as is defined by the alphabet.
     * E.g. for DNA, this will map {0 -> 'A', 1 -> 'C', ...}.
     * This function also works
     * in-place, i.e. `first` can be equal to `result`. In the case of overlapping
     * sequences with `first` != `result`, the behavious is not defined.
     *
     * @param first   InputIterator to the start of the sequence to be translated.
     * @param last    InputIterator to the end of the sequence to be translated
     *                (that is one element past the last element in the sequence).
     * @param result  OuptutIterator to the start of the result sequence.
     *
     * @tparam InputIterator  InputIterator type for the `begin()` and `end()`
     *                        iterators.
     * @tparam OutputIterator Iterator type for sequential output operations.
     *
     * @return        An Iterator to the end of the result sequence, this will
     *                point to one element past the last written element.
     */
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

  } // namespace common
} // namespace bliss

#endif // BLISS_COMMON_ALPHABETTRAITS_H
