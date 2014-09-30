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
#ifndef BLISS_COMMON_PACKINGITERATORS_H
#define BLISS_COMMON_PACKINGITERATORS_H

// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <iterator>
#include <type_traits>

// own includes
#include <common/base_types.hpp>
#include <common/bit_ops.hpp>
#include <common/padding.hpp>
#include <common/Kmer.hpp>

#include <iterators/many2one_iterator.hpp>
#include <iterators/one2many_iterator.hpp>

namespace bliss
{

/**
 * @brief Functor for the char packing iterator.
 *
 * Packs characters of a certain bit size into the given storage type.
 * E.g. packing 16 characters of 2 bits each into a 32bit word.
 *
 * @tparam BITS_PER_CHAR
 * @tparam PackedStorageType
 */
template <BitSizeType BITS_PER_CHAR, typename PackedStorageType=WordType>
struct PackingItFunctor
{
  /// the properties of potential padding in the storage words, given
  /// the storage type and the bits per character.
  typedef PackingTraits<PackedStorageType, BITS_PER_CHAR> padtraits;
  /// The number of characters packed into each word
  static constexpr unsigned int m = padtraits::chars_per_word;

  /**
   * @brief Packs characters from the given iterator `cur` into a single
   *        storage word and returns that.
   *
   * Reads as many characters from the `cur` iterator as are needed to fill a
   * single storage word of type `PackedStorageType`. This corresponds to at
   * most `m` characters. If the `cur` iterator reaches `end`, no more
   * characters are read and the return value is correspondingly padded with
   * `0` characters.
   *
   * @tparam Iterator   The type of the iterator returning characters.
   * @param cur[in|out] The current iterator position. The needed number of
   *                    characters are read from this iterator and the iterator
   *                    is set to the position of the character that needs to
   *                    be read next.
   * @param end         Stop reading characters from the input iterator as soon
   *                    as the input iterator `cur` compares equal to this
   *                    `end` iterator.
   *
   * @return The packed storage word containing all read characters.
   */
  template<typename Iterator>
  PackedStorageType operator()(Iterator& cur, const Iterator& end)
  {
    // get underlying type of iterator
    typedef typename std::iterator_traits<Iterator>::value_type base_value_type;
    // advance the next iterator and buffer result
    BitSizeType numChars = static_cast<BitSizeType>(std::min<std::size_t>(m, std::distance(cur, end)));
    std::array<base_value_type, m> word_cache;
    for (BitSizeType i = 0; i < numChars; ++i)
    {
      word_cache[i] = *cur;
      ++cur;
    }
    // read them in reverse
    PackedStorageType buffer = 0;
    for (BitSizeType i = numChars; i > 0; --i)
      //for (BitSizeType i = 1; i <= numChars; ++i)
    {
      buffer <<= padtraits::bits_per_char;
      buffer |= getLeastSignificantBitsMask<PackedStorageType>(padtraits::bits_per_char) & word_cache[i-1];
    }

    return buffer;
  }
};

/**
 * @brief Packs multiple characters from an underlying base iterator into
 *        larger words.
 *
 * @tparam BaseIterator         The type of the underlying iterator.
 * @tparam bits_per_char        The number of data bits per character.
 * @tparam PackedStorageType    The storage type.
 */
template <typename BaseIterator, BitSizeType bits_per_char, typename PackedStorageType=WordType>
class PackingIterator
  : public iterator::many2one_iterator<BaseIterator,
                      PackingItFunctor<bits_per_char, PackedStorageType> >
{
  /// The base class type, i.e. the many2one iterator templated to this
  /// use-case
  typedef iterator::many2one_iterator<BaseIterator,
            PackingItFunctor<bits_per_char, PackedStorageType> > base_class_t;
  /// The operator type PackingItFunctor, templated to the parameters of this
  /// class
  typedef PackingItFunctor<bits_per_char, PackedStorageType> functor_t;
  /// The padding traits (number of padding bits and data bits, etc)
  typedef PackingTraits<PackedStorageType, bits_per_char> padtraits;

public:
  /// Default constructor
  PackingIterator() {}

  /**
   * @brief Constructs a non-functional iterator of the given base position.
   *
   * This will set the base iterator `begin` and `end` position to the given
   * iterator position. Thus, this iterator can not be dereferenced or
   * advanced any further. This can be used to construct an `end` iterator
   * that is solely used for comparison against a functional iterator.
   *
   * @param end   The iterator position to fix this iterator to.
   */
  PackingIterator(BaseIterator end)
    : base_class_t(end, end, functor_t(), padtraits::chars_per_word)  {}

  /**
   * @brief Constructs a packing iterator based on the given base iterator.
   *
   * @param baseBegin   The iterator position to begin reading characters from.
   * @param baseEnd     The `end` iterator for the given input iterator range.
   *                    Characters will be read from the baseBegin iterator
   *                    until it compares equal to this iterator.
   */
  PackingIterator(BaseIterator baseBegin, BaseIterator baseEnd)
    : base_class_t(baseBegin, baseEnd, functor_t(), padtraits::chars_per_word)  {}
};

/**
 * @brief Functor for the unpacking iterator.
 *
 * Unpacks a given word into a sequence of characters.
 *
 * @tparam bits_per_char    The number of bits per character.
 * @tparam unpacked_type    The unpacked type (most likely = `char`)
 */
template <unsigned int bits_per_char, typename unpacked_type=CharType>
struct UnpackingFunctor
{
  /**
   * @brief Unpacks one character at a time from the given word `value` from
   *        at the given offset `offset`.
   *
   * @tparam T          The value (= packed word) type.
   * @tparam diff_type  The offset type.
   * @param value       The storage word to unpack.
   * @param offset      The offset inside the word from where to unpack the
   *                    next character.
   * @return The unpacked character.
   */
  template <typename T, typename diff_type>
  unpacked_type operator()(const T& value, const diff_type& offset)
  {
    // directly access the according position, no need to buffer
    T mask = getLeastSignificantBitsMask<T>(bits_per_char);
    // shift and mask to get the value
    T x = mask & (value >> (offset*bits_per_char));
    // cast to target type
    return static_cast<unpacked_type>(x);
  }
};

/**
 * @brief Unpacks characters from storage words given by an underlying iterator.
 *
 * @tparam BaseIterator         The type of the underlying iterator.
 * @tparam bits_per_char        The number of data bits per character.
 * @tparam UnpackedStorageType  The data type of the unpacked characters.
 */
template <typename BaseIterator, BitSizeType bits_per_char, typename UnpackedStorageType=CharType>
class UnpackingIterator
  : public iterator::one2many_iterator<BaseIterator,
            UnpackingFunctor<bits_per_char, UnpackedStorageType> >
{
  /// The iterator base type (one2many) templated to this use case
  typedef iterator::one2many_iterator<BaseIterator,
            UnpackingFunctor<bits_per_char, UnpackedStorageType> > base_class_t;
  /// The value type of the underlying base iterator
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;

  /// The difference type of this iterator
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;
  /// The type of the functor used for unpacking
  typedef UnpackingFunctor<bits_per_char, UnpackedStorageType> functor_t;
  /// The padding properties (data bits, padding bits, etc)
  typedef PackingTraits<base_value_type, bits_per_char> padtraits;

public:
  /// Default constructor
  UnpackingIterator() {}

  /**
   * @brief Constructs an unpacking iterator starting from the given base
   *        iterator position.
   *
   * @param baseBegin The base iterator position to start unpacking from.
   */
  UnpackingIterator(BaseIterator baseBegin)
    : base_class_t(baseBegin, functor_t(), padtraits::chars_per_word)  {}

  /**
   * @brief Constructs an unpacking iterator starting from the given base
   *        iterator position advanced by `offset` number of characters.
   *
   * This constructor advances the iterator position by `offset` many
   * characters. This can be used to create `end` style iterator, when the
   * number of characters to read from the base iterator is known in advance.
   * NOTE that this is necessary if the last storage word in the input sequence
   * contains less than the full amount of packed characters, since there is
   * no encoding of how many characters are packed in each storage word.
   *
   * @param baseBegin   The base iterator position to start from
   * @param offset      An additional offset in number of characters which are
   *                    skipped in the input iterator during construction.
   */
  UnpackingIterator(BaseIterator baseBegin, diff_type offset)
    : base_class_t(baseBegin, functor_t(), padtraits::chars_per_word, offset)  {}
};

} // namespace bliss

#endif // BLISS_COMMON_PACKINGITERATORS_H
