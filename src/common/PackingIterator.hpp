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
#include <type_traits>

// own includes
#include <common/base_types.hpp>
#include <common/bit_ops.hpp>
#include <common/padding.hpp>
#include <common/Kmer.hpp>

#include <iterators/many2one_iterator.hpp>
#include <iterators/one2many_iterator.hpp>
#include <iterators/sliding_window_iterator.hpp>

namespace bliss
{

/*

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
*/


template <BitSizeType BITS_PER_CHAR, typename PackedStorageType=WordType>
struct PackingItFunctor
{

  typedef PaddingTraits<PackedStorageType, BITS_PER_CHAR> padtraits;
  static constexpr unsigned int m = padtraits::chars_per_word;

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
        buffer |= getBitMask<PackedStorageType>(padtraits::bits_per_char) & word_cache[i-1];
      }

      return buffer;
    }
};

/**
 * @brief
 */
template <typename BaseIterator, BitSizeType bits_per_char, typename PackedStorageType=WordType>
class PackingIterator : public iterator::many2one_iterator<BaseIterator, PackingItFunctor<bits_per_char, PackedStorageType> >
{
  typedef iterator::many2one_iterator<BaseIterator, PackingItFunctor<bits_per_char, PackedStorageType> > base_class_t;
  typedef PackingItFunctor<bits_per_char, PackedStorageType> functor_t;
  typedef PaddingTraits<PackedStorageType, bits_per_char> padtraits;

public:
  PackingIterator() {}

  PackingIterator(BaseIterator baseBegin)
    : base_class_t(baseBegin, baseBegin, functor_t(), padtraits::chars_per_word)  {}

  PackingIterator(BaseIterator baseBegin, BaseIterator baseEnd)
    : base_class_t(baseBegin, baseEnd, functor_t(), padtraits::chars_per_word)  {}
};

template <unsigned int bits_per_char, typename unpacked_type=CharType>
struct UnpackingFunctor
{
  template <typename T, typename diff_type>
  unpacked_type operator()(const T& value, const diff_type& offset)
  {
    // directly access the according position, no need to buffer
    T mask = getBitMask<T>(bits_per_char);
    // shift and mask to get the value
    T x = mask & (value >> (offset*bits_per_char));
    // cast to target type
    return static_cast<unpacked_type>(x);
  }
};

template <typename BaseIterator, BitSizeType bits_per_char, typename UnpackedStorageType=CharType>
class UnpackingIterator : public iterator::one2many_iterator<BaseIterator, UnpackingFunctor<bits_per_char, UnpackedStorageType> >
{
  typedef iterator::one2many_iterator<BaseIterator, UnpackingFunctor<bits_per_char, UnpackedStorageType> > base_class_t;
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;
  typedef UnpackingFunctor<bits_per_char, UnpackedStorageType> functor_t;
  typedef PaddingTraits<base_value_type, bits_per_char> padtraits;

public:
  UnpackingIterator() {}

  UnpackingIterator(BaseIterator baseBegin)
    : base_class_t(baseBegin, functor_t(), padtraits::chars_per_word)  {}

  UnpackingIterator(BaseIterator baseBegin, diff_type offset)
    : base_class_t(baseBegin, functor_t(), padtraits::chars_per_word, offset)  {}
};

/*
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
    // directly access the according position, no need to buffer
    base_value_type mask = getBitMask<base_value_type>(bitsPerChar);
    base_value_type x = mask & (*baseCur >> (curCharOffset*bitsPerChar));
    // cast to target type
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
*/

/**
 * @brief The sliding window operator for k-mer generation from character data.
 *
 * @tparam BaseIterator Type of the underlying base iterator, which returns
 *                      characters.
 * @tparam Kmer         The k-mer type, must be of type bliss::Kmer
 */
template <class BaseIterator, class Kmer>
class KmerSlidingWindow {};

template <typename BaseIterator, unsigned int KMER_SIZE,
          unsigned int BITS_PER_CHAR, typename word_type>
class KmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> >
{
public:
  /// The Kmer type (same as the `value_type` of this iterator)
  typedef bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> kmer_type;
  /// The value_type of the underlying iterator
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;

  /**
   * @brief Initializes the sliding window.
   *
   * @param it[in|out]  The current base iterator position. This will be set to
   *                    the last read position.
   */
  inline void init(BaseIterator& it)
  {
    kmer.fillFromChars(it, true);
  }

  /**
   * @brief Slides the window by one character taken from the given iterator.
   *
   * This will read the current character of the iterator and then advance the
   * iterator by one.
   *
   * @param it[in|out]  The underlying iterator position, this will be read
   *                    and then advanced.
   */
  inline void next(BaseIterator& it)
  {
    kmer.nextFromChar(*it);
    ++it;
  }

  /**
   * @brief Returns the value of the current sliding window, i.e., the current
   *        k-mer value.
   *
   * @return The current k-mer value.
   */
  inline kmer_type getValue()
  {
    // return a copy of the current kmer
    return this->kmer;
  }
private:
  /// The kmer buffer (i.e. the window of the sliding window)
  kmer_type kmer;
};


/**
 * @brief Iterator that generates k-mers from character data.
 *
 * @tparam BaseIterator     The underlying iterator of characters. E.g. a
 *                          std::string::iterator. Any iterator yielding
 *                          `char` works.
 * @tparam Kmer             The type of the Kmer, this has to be of type
 *                          bliss::Kmer.
 */
template <class BaseIterator, class Kmer>
class KmerGenerationIterator {};

// template specialization for bliss::Kmer as kmer type
// The template parameters KMER_SIZE, BITS_PER_CHAR and word_type are
// set to the appropriate template parameters of bliss::Kmer and do not
// have to be explicitly stated when creating this class.
template <typename BaseIterator, unsigned int KMER_SIZE, unsigned int BITS_PER_CHAR, typename word_type>
class KmerGenerationIterator<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> >
: public iterator::sliding_window_iterator<BaseIterator, KmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> > >
{
protected:
  /// The type of the base class
  typedef iterator::sliding_window_iterator<BaseIterator, KmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> > >
 base_class_t;

  /// The difference_type of character offsets
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;

  /// The type of the sliding window.
  typedef KmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> > functor_t;

public:
  /// Default constructor.
  KmerGenerationIterator() : base_class_t() {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first character of the
   *                      sequence to be used for generating k-mers.
   */
  KmerGenerationIterator(const BaseIterator& baseBegin)
    : base_class_t(baseBegin, true) {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param baseBegin     An iterator pointing to the first character of the
   *                      sequence to be used for generating k-mers.
   * @param initialze_window    Whether to read the first `KMER_SIZE` characters
   *                            to initialize the sliding window. This will
   *                            have to be set to `false` for creating
   *                            `end` style iterators which would otherwise
   *                            read past the valid range.
   */
  KmerGenerationIterator(const BaseIterator& baseBegin, bool initialize_window)
    : base_class_t(baseBegin, initialize_window) {}

protected:
  /*****************************
   *  non public constructors  *
   *****************************/
  // handleing of the `window` object instances is strictly hidden

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param baseBegin     An iterator pointing to the first character of the
   *                      sequence to be used for generating k-mers.
   * @param window        The sliding window object used by the iterator.
   */
  KmerGenerationIterator(const BaseIterator& baseBegin, const functor_t& window)
    : base_class_t(baseBegin, window, true)  {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first character of the
   *                      sequence to be used for generating k-mers.
   * @param window        The sliding window object used by the iterator.
   * @param initialze_window    Whether to read the first `KMER_SIZE` characters
   *                            to initialize the sliding window. This will
   *                            have to be set to `false` for creating
   *                            `end` style iterators which would otherwise
   *                            read past the valid range.
   */
  KmerGenerationIterator(const BaseIterator& baseBegin, const functor_t& window, bool initialize_window)
    : base_class_t(baseBegin, window, initialize_window)  {}
};


/**
 * @brief The sliding window operator for k-mer generation from packed data.
 *
 * @tparam BaseIterator Type of the underlying base iterator, which returns
 *                      elements of packed data.
 * @tparam Kmer         The k-mer type, must be of type bliss::Kmer
 */
template <class BaseIterator, class Kmer>
class PackedKmerSlidingWindow {};

template <typename BaseIterator, unsigned int KMER_SIZE,
          unsigned int BITS_PER_CHAR, typename word_type>
class PackedKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> >
{
public:
  /// The Kmer type (same as the `value_type` of this iterator)
  typedef bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> kmer_type;
  /// The value_type of the underlying iterator
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;
  /// The padding traits of the underlying stream
  typedef PaddingTraits<base_value_type, BITS_PER_CHAR> padtraits;

  /**
   * @brief Initializes the sliding window.
   *
   * @tparam offset_t       The offset type.
   * @param it[in|out]      The current base iterator, this will be set
   *                        to the last read position.
   * @param offset[in|out]  The current bit/character offset in the current position
   *                        of the base iterator. This will be set to the
   *                        position that was read last.
   */
  template<typename offset_t>
  inline void init(BaseIterator& it, offset_t& offset)
  {
    // there is no implementation yet to handle the case that the first
    // k-mer starts from an offset != 0:
    assert(offset == 0);

    // fill kmer from the given packed and padded stream
    // this leaves the iterator and the offset ON the last read position
    // and NOT AFTER this (i.e. NOT on the position to be read NEXT)
    offset = kmer.fillFromPaddedStream(it, true);
  }

  /**
   * @brief Slides the window by one character, taken from the given iterator.
   *
   * This will read the current character of the iterator and then increase
   * the (iterator, offset) position by one character.
   *
   * @tparam offset_t       The offset type.
   * @param it[in|out]      The current base iterator, this will be set
   *                        to the next position.
   * @param offset[in|out]  The current bit/character offset in the current position
   *                        of the base iterator. This will be set to one
   *                        character past the character read.
   */
  template<typename offset_t>
  inline void next(BaseIterator& it, offset_t& offset)
  {
    kmer.nextFromPaddedStream(it, offset);
  }

  /**
   * @brief Returns the value of the current sliding window, i.e., the current
   *        k-mer value.
   *
   * @return The current k-mer value.
   */
  inline kmer_type getValue()
  {
    // return a copy of the current kmer
    return this->kmer;
  }

  /**
   * @brief Skips over `advance_by` characters in the current (iterator,
   *        offset) position.
   *
   * @tparam offset_t   Type of the offset.
   * @param it[in|out]      The current iterator position. Will be modified to
   *                        the position after skipping `advance_by` chars.
   * @param offset          The current offset position. Will be modified to
   *                        the position after skipping `advance_by` chars.
   * @param advance_by      The number of characters to skip.
   */
  template<typename offset_t>
  // TODO: - [ ] separate offset types
  //         [ ] FIX: from any starting offset, not only from offset=0
  //
  void skip(BaseIterator& it, offset_t& offset, offset_t advance_by)
  {
    // get offset for underlying iterator
    offset_t nWords = advance_by / padtraits::chars_per_word;
    // get offset for bitwise offset
    offset_t bit_offset = advance_by % padtraits::chars_per_word;
    bit_offset *= padtraits::bits_per_char;

    // add it to the base iterator
    std::advance(it, nWords);
    offset += bit_offset;
  }
private:
  /// The kmer buffer (i.e. the window of the sliding window)
  kmer_type kmer;
};



/**
 * @brief Iterator that generates k-mers from packed data.
 *
 * @tparam BaseIterator     The underlying iterator, which supplies any integer
 *                          type of packed character data.
 * @tparam Kmer             The type of the Kmer, this has to be of type
 *                          bliss::Kmer.
 */
template <class BaseIterator, class Kmer>
class PackedKmerGenerationIterator {};

// template specialization for bliss::Kmer as kmer type
// The template parameters KMER_SIZE, BITS_PER_CHAR and word_type are
// set to the appropriate template parameters of bliss::Kmer and do not
// have to be explicitly stated when creating this class.
template <typename BaseIterator, unsigned int KMER_SIZE,
          unsigned int BITS_PER_CHAR, typename word_type>
class PackedKmerGenerationIterator<BaseIterator,
      bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> >
: public iterator::one2many_sliding_window_iterator<BaseIterator,
          PackedKmerSlidingWindow<BaseIterator,
                bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> > >
{
protected:
  /// The type of the base class
  typedef iterator::one2many_sliding_window_iterator<BaseIterator, PackedKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> > >
 base_class_t;

  /// The difference_type of character offsets
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;

  /// The type of the sliding window.
  typedef PackedKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> > functor_t;

public:
  /// Default constructor
  PackedKmerGenerationIterator() : base_class_t() {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first element of the
   *                      packed and padded sequence to be used for generating
   *                      k-mers.
   */
  PackedKmerGenerationIterator(const BaseIterator& baseBegin)
    : base_class_t(baseBegin, functor_t())  {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point and advancing
   *          `char_offset` characters into the underlying sequence.
   *
   * @param   baseBegin   An iterator pointing to the first element of the
   *                      packed and padded sequence to be used for generating
   *                      k-mers.
   * @param  offset       A character offset, the constructed iterator is
   *                      advanced by this many positions at construction time.
   *                      This can be used to create the `end` iterator for
   *                      generating sequences.
   */
  PackedKmerGenerationIterator(const BaseIterator& baseBegin, diff_type offset)
    : base_class_t(baseBegin, functor_t(), offset)  {}

protected:
  /*****************************
   *  non public constructors  *
   *****************************/
  // handleing of the `window` object instances is strictly hidden

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first element of the
   *                      packed and padded sequence to be used for generating
   *                      k-mers.
   * @param window        The sliding window object used by the iterator.
   */
  PackedKmerGenerationIterator(const BaseIterator& baseBegin,
                               const functor_t& window)
    : base_class_t(baseBegin, window)  {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point and advancing
   *          `char_offset` characters into the underlying sequence.
   *
   * @param   baseBegin   An iterator pointing to the first element of the
   *                      packed and padded sequence to be used for generating
   *                      k-mers.
   * @param  offset       A character offset, the constructed iterator is
   *                      advanced by this many positions at construction time.
   *                      This can be used to create the `end` iterator for
   *                      generating sequences.
   * @param window        The sliding window object used by the iterator.
   */
  PackedKmerGenerationIterator(const BaseIterator& baseBegin,
                               const functor_t& window, diff_type offset)
    : base_class_t(baseBegin, window, offset)  {}
};

#if 0
// TODO add template specialization for padding (TRUE or FALSE)
template <typename BaseIterator, unsigned int KMER_SIZE, unsigned int BITS_PER_CHAR, typename word_type>
class KmerGenerationIterator<BaseIterator, bliss::Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> >
  : public std::iterator<std::forward_iterator_tag,//typename std::iterator_traits<BaseIterator>::iterator_category,
                  Kmer<KMER_SIZE, BITS_PER_CHAR, word_type>,
                  typename std::iterator_traits<BaseIterator>::difference_type>
{
public:
  /// The Kmer type (same as the `value_type` of this iterator)
  typedef Kmer<KMER_SIZE, BITS_PER_CHAR, word_type> kmer_type;
  /// The value type of the iterator. This is the Kmer type.
  typedef kmer_type value_type;

  /// The type of the underlying base iterator.
  typedef BaseIterator base_iterator;
  /// The value type of the underlying base iterator.
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;

  /// The bits per character in the k-mer representation
  static constexpr BitSizeType bitsPerChar = BITS_PER_CHAR;

  /**
   * @brief   Constructs a new KmerGenerationIterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first element of the
   *                      packed and padded sequence to be used for generating
   *                      k-mers.
   */
  KmerGenerationIterator(BaseIterator baseBegin)
    : baseCurrent(baseBegin), bitOffsetCurrent(0),
      baseRead(baseBegin), bitOffsetRead(0)
  {
    initFirstKmer();
  }
  /**
   * @brief   Constructs a new KmerGenerationIterator using the underlying
   *          iterator `baseBegin` as starting point and advancing
   *          `char_offset` characters into the underlying sequence.
   *
   * @param   baseBegin   An iterator pointing to the first element of the
   *                      packed and padded sequence to be used for generating
   *                      k-mers.
   * @param   char_offset A character offset, this iterator is advanced by
   *                      this many positions at construction time. This
   *                      can be used to create the `end` iterator for
   *                      generating sequences.
   */
  KmerGenerationIterator(BaseIterator baseBegin, unsigned int char_offset)
    : baseCurrent(baseBegin), bitOffsetCurrent(0),
      baseRead(baseBegin), bitOffsetRead(0)
  {
    advance_chars(char_offset);

    // TODO: initFirstKmer with offset !?
    //      -> right now this is only used as END iterator that is never increased
    //         or dereferenced
  }

  /**
   * @brief   Returns the current k-mer.
   *
   * @returns A copy of the current k-mer.
   */
  value_type operator*()
  {
    while (incOffset > 0)
    {
      // read through sliding window
      // TODO: (enhancement) for big `incOffset` it should be more efficient to
      //       construct the kmer from scratch
      kmer_buffer.nextFromPaddedStream(baseRead, bitOffsetRead);
      --incOffset;
    }

    // return a copy of the kmer
    return kmer_buffer;
  }

  /**
   * @brief   Iterates this iterator to one higher position.
   *
   * @returns A reference to this object.
   */
  KmerGenerationIterator& operator++() // prefix ++
  {
    // increment the comparable base iterator
    bitOffsetCurrent += bitsPerChar;
    if (bitOffsetCurrent == PaddingTraits<base_value_type, bitsPerChar>::data_bits)
    {
      ++baseCurrent;
      bitOffsetCurrent = 0;
    }

    // increase the read offset
    ++incOffset;

    // return the obligatory self reference
    return (*this);
  }

  /**
   * @brief   Returns whether this and the given iterator are equal.
   *
   * Two iterators compare equal if the underlying iterators compare equal and
   * the current offset is identical.
   *
   * @param rhs   The iterator to compare this iterator against.
   * @returns     `true` if the iterators compare equal, `false` otherwise.
   */
  bool operator==(const KmerGenerationIterator& rhs)
  {
    // compare the base iterator and its associated offset
    return (this->baseCurrent == rhs.baseCurrent) && (this->bitOffsetCurrent == rhs.bitOffsetCurrent);
  }

  /**
   * @brief   Returns whether this and the given iterator are not equal.
   *
   * Two iterators compare equal if the underlying iterators compare equal and
   * the current offset is identical. Otherwise they are not equal, and this
   * function will return `true`.
   *
   * @param rhs   The iterator to compare this iterator against.
   * @returns     `true` if the iterators are not equal, `false` otherwise.
   */
  bool operator!=(const KmerGenerationIterator& rhs)
  {
    return !this->operator==(rhs);
  }

  /**
   * @brief The postfix increase operator.
   * @see   operator++()
   *
   */
  KmerGenerationIterator& operator++(int) // postfix ++
  {
    // create a copy prior to increasing
    KmerGenerationIterator tmp(*this);
    this->operator++();
    return tmp;
  }

private:

  /**
   * @brief   Initializes the first k-mer from the base iterator and
   *          advances the corresponding iterators.
   */
  void initFirstKmer()
  {
    bitOffsetRead = kmer_buffer.fillFromPaddedStream(baseRead, true);
    // init the `Current` iterator at one position prior
    // TODO: replace KMER_SIZE by Kmer.getSize() operation
    advance_chars(KMER_SIZE);
    // dumb workaround
    std::swap(bitOffsetRead, bitOffsetCurrent);
    std::swap(baseRead, baseCurrent);
  }

  /**
   * @brief   Advances the `Current` iterator and its associated bit offset
   *          by `char_offset` many characters.
   *
   * @param char_offset   The character offset to advance by.
   */
  void advance_chars(unsigned int char_offset)
  {
    // get the word and bit offset
    unsigned int nWords = char_offset / PaddingTraits<base_value_type, bitsPerChar>::chars_per_word;
    unsigned int bit_offset = char_offset % PaddingTraits<base_value_type, bitsPerChar>::chars_per_word;
    bit_offset *= bitsPerChar;

    // add it to the base iterator
    std::advance(baseCurrent, nWords);
    bitOffsetCurrent += bit_offset;
  }


  /*
   * This iterator holds two independent copies of the base iterator:
   *
   * The `Current` iterator and its according offset mirror the current
   * position of the iterator. This is only used to compare two copies
   * of this type to one another, and is especially NOT used to read the
   * underlying data.
   *
   * The `Read` iterator and its according offset are used as iterators
   * for the sliding window k-mer generation. This is used to read the data
   * from the base iterator.
   */

  // The `Current` iterator and its associated bit offset
  /// The underlying base iterator pointing to the current position.
  BaseIterator baseCurrent;
  /// The bit offset in the current word pointed to by `baseCurrent`.
  unsigned int bitOffsetCurrent = 0;

  // The `Read` iterator and its associated bit offset
  /// The underlying base iterator pointing to the next position to be read.
  BaseIterator baseRead;
  /// The bit offset in the word pointed to by the previous iterator.
  unsigned int bitOffsetRead = 0;

  // The number of increases that the `Read` iterator laggs behind
  // NOTE: this is the number of window shifts have to be done on the k-mer
  //       and NOT the number if increases that have to be done to compare
  //       equal to the `Current` iterator
  /// The lagg between the `current` and `read` iterators.
  unsigned int incOffset = 0;

  /// The current kmer
  kmer_type kmer_buffer;
};

#endif

} // namespace bliss

#endif // BLISS_COMMON_PACKINGITERATOR_H
