/**
 * @file    kmer_iterators.hpp
 * @ingroup common
 * @author  Patrick Flick
 * @brief
 *
 * Copyright (c) TODO
 *
 * TODO add Licence
 */
#ifndef BLISS_COMMON_KMER_ITERATORS_H
#define BLISS_COMMON_KMER_ITERATORS_H

// C std lib includes:
#include <cstdlib>

// C++ STL includes:
#include <iterator>
#include <type_traits>

// own includes
#include <common/base_types.hpp>
#include <common/padding.hpp>
#include <common/Kmer.hpp>

#include <iterators/sliding_window_iterator.hpp>

// TODO: Need convenience functions to make start and end iterators (the true/false flags are not good.  else enforce that flag has no default value.)
// TODO: Need convenience typedef to allow for reverse complement generation.

namespace bliss
{

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
          typename ALPHABET, typename word_type>
class KmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, ALPHABET, word_type> >
{
public:
  /// The Kmer type (same as the `value_type` of this iterator)
  typedef bliss::Kmer<KMER_SIZE, ALPHABET, word_type> kmer_type;
  typedef BaseIterator  base_iterator_type;
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
 * @brief The sliding window operator for reverse k-mer generation from character data.
 * @note  to create reverse complement, use a transform iterator to change the char to complement
 *
 * @tparam BaseIterator Type of the underlying base iterator, which returns
 *                      characters.
 * @tparam Kmer         The k-mer type, must be of type bliss::Kmer
 */
template <class BaseIterator, class Kmer>
class ReverseKmerSlidingWindow {};

template <typename BaseIterator, unsigned int KMER_SIZE,
          typename ALPHABET, typename word_type>
class ReverseKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, ALPHABET, word_type> >
{
public:
  /// The Kmer type (same as the `value_type` of this iterator)
  typedef bliss::Kmer<KMER_SIZE, ALPHABET, word_type> kmer_type;
  typedef BaseIterator  base_iterator_type;
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
    kmer.fillReverseFromChars(it, true);
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
    kmer.nextReverseFromChar(*it);
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
template <class SlidingWindow>
class KmerGenerationIteratorBase
    : public iterator::sliding_window_iterator<typename SlidingWindow::base_iterator_type, SlidingWindow >
{
protected:
    typedef typename SlidingWindow::base_iterator_type BaseIterator;

  /// The type of the base class
  typedef iterator::sliding_window_iterator<BaseIterator, SlidingWindow >
 base_class_t;

  /// The difference_type of character offsets
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;

  /// The type of the sliding window.
  typedef SlidingWindow functor_t;

public:
  /// Default constructor.
  KmerGenerationIteratorBase() : base_class_t() {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first character of the
   *                      sequence to be used for generating k-mers.
   */
  KmerGenerationIteratorBase(const BaseIterator& baseBegin)
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
  KmerGenerationIteratorBase(const BaseIterator& baseBegin, bool initialize_window)
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
  KmerGenerationIteratorBase(const BaseIterator& baseBegin, const functor_t& window)
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
  KmerGenerationIteratorBase(const BaseIterator& baseBegin, const functor_t& window, bool initialize_window)
    : base_class_t(baseBegin, window, initialize_window)  {}
};

template <class BaseIterator, class Kmer>
using KmerGenerationIterator = KmerGenerationIteratorBase<KmerSlidingWindow<BaseIterator, Kmer > >;

template <class BaseIterator, class Kmer>
using ReverseKmerGenerationIterator = KmerGenerationIteratorBase<ReverseKmerSlidingWindow<BaseIterator, Kmer > >;




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
          typename ALPHABET, typename word_type>
class PackedKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, ALPHABET, word_type> >
{
public:
  /// The Kmer type (same as the `value_type` of this iterator)
  typedef bliss::Kmer<KMER_SIZE, ALPHABET, word_type> kmer_type;
  /// The value_type of the underlying iterator
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;
  /// The padding traits of the underlying stream
  typedef PackingTraits<base_value_type, bliss::AlphabetTraits<ALPHABET>::getBitsPerChar()> padtraits;

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
    offset = kmer.fillFromPackedStream(it, offset, true);
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
    kmer.nextFromPackedStream(it, offset);
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
// The template parameters KMER_SIZE, ALPHABET and word_type are
// set to the appropriate template parameters of bliss::Kmer and do not
// have to be explicitly stated when creating this class.
template <typename BaseIterator, unsigned int KMER_SIZE,
          typename ALPHABET, typename word_type>
class PackedKmerGenerationIterator<BaseIterator,
      bliss::Kmer<KMER_SIZE, ALPHABET, word_type> >
: public iterator::one2many_sliding_window_iterator<BaseIterator,
          PackedKmerSlidingWindow<BaseIterator,
                bliss::Kmer<KMER_SIZE, ALPHABET, word_type> > >
{
protected:
  /// The type of the base class
  typedef iterator::one2many_sliding_window_iterator<BaseIterator, PackedKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, ALPHABET, word_type> > >
 base_class_t;

  /// The difference_type of character offsets
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;

  /// The type of the sliding window.
  typedef PackedKmerSlidingWindow<BaseIterator, bliss::Kmer<KMER_SIZE, ALPHABET, word_type> > functor_t;

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

} // namespace bliss

#endif // BLISS_COMMON_KMER_ITERATORS_H
