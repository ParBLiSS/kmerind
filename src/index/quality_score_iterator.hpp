/**
 * @file    quality_score_iterator.hpp
 * @ingroup index
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @author  Tony Pan <tcp1975@gmail.com>
 *
 * @brief   Implements the sliding window quality score generating iterator.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 *
 * TODO add Licence
 */

#ifndef BLISS_INDEX_QUALITY_SCORE_ITERATOR_HPP
#define BLISS_INDEX_QUALITY_SCORE_ITERATOR_HPP

#include <vector>

#include <index/quality_scores.hpp>
#include <iterators/sliding_window_iterator.hpp>

namespace bliss
{
namespace index
{

/**
 * @brief The sliding window operator for calculating quality scores for k-mers.
 *
 * @details
 *
 * If `p_i` is the probability of a base `i` to be correct, then the probability
 * of a k-mer with bases `b_1,...,b_k` being correct is given as
 *     `p_1 * p_2 * ... * p_k`
 *
 * Since multiplications of floating point numbers is inefficient, we calculate
 * the product in log_2 space in the sliding window. Thus we simply sum the
 * log probabilities into a common sum. Sliding the window by one position thus
 * means that we have to subtract the value "falling out of the window" and
 * add the new value coming in.
 *
 * Retrieving the k-mer probability is then done by calculating the exponent
 * of the current sum to the base of 2.
 *
 * @tparam BaseIterator Type of the underlying base iterator, which returns
 *                      quality scores from the reads.
 * @tparam Encoder     The `Encoder` class, which decodes/encodes
 *                      between specific quality score (in ASCII space) and
 *                      the actual probabilities of bases being correct
 *                      (probabilities in log space).
 */
template <typename BaseIterator, unsigned int KMER_SIZE,
          typename Encoder = bliss::index::Illumina18QualityScoreCodec<float> >
class QualityScoreSlidingWindow
{
private:
  /// Type of the internal sum is equal to the floating point type
  /// returned by the encoding
  typedef typename Encoder::value_type QualityType;

  /// The current sum of base-correct probabilities in the current window
  QualityType current_sum = 0;

  /// All current values in the window, updated as a circular buffer
  /// since the order does not matter.
  QualityType window_values[KMER_SIZE];

  /// number of incorrect bases in the current window (i.e., bases with zero
  /// probability of being correct)
  unsigned int n_incorrect_bases = 0;

  /// The current position in the window_values circular buffer
  /// Points to the next available position
  unsigned int window_pos = 0;
public:
  /// The value_type of the underlying iterator
  typedef typename std::iterator_traits<BaseIterator>::value_type base_value_type;

  // TODO: constructors

  /**
   * @brief Initializes the sliding window.  at end of initialization, it is set to last position read.
   *
   * @param it[in|out]  The current base iterator position. This will be set to
   *                    the last read position (and NOT the next one)
   */
  inline void init(BaseIterator& it)
  {
    for (unsigned int i = 0; i < KMER_SIZE;)
    {
      // fill the window_values and count the number of bases with
      // zero probability (-inf log-prob of) of begin correct
      QualityType newval = Encoder::decode(*it);
      window_values[i] = newval;

      if (newval < Encoder::min_log_prob)
      {
        ++n_incorrect_bases;
      } else {
        current_sum += newval;
      }
      if (++i < KMER_SIZE) ++it;
    }
    window_pos = 0;
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
    // slide the window by one, throwing out the last value and adding in
    // the new one in a circular buffer fashion
    QualityType oldval = window_values[window_pos];
    QualityType newval = Encoder::decode(*it);

    window_values[window_pos] = newval;
    window_pos = (window_pos+1) % KMER_SIZE;

    // remove old value from either the sum or the incorrect count
    if (oldval < Encoder::min_log_prob)
    {
      --n_incorrect_bases;
    }
    else
    {
      current_sum -= oldval;
    }

    // add the new value to either the sum or the incorrect count
    if (newval < Encoder::min_log_prob)
    {
      ++n_incorrect_bases;
    }
    else
    {
      current_sum += newval;
    }
    ++it;
  }

  /**
   * @brief Returns the probability [0,1] that the k-mer contains only correct
   *        bases. This is equal to the multiplication of all base probabilities
   *        of being correct.
   *
   * @return The probability that all bases in the sliding window are correct.
   */
  inline QualityType getValue()
  {

      if (n_incorrect_bases > 0)
        return 0.0;
      else
        return std::exp2(current_sum);
  }
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
template <typename BaseIterator, unsigned int KMER_SIZE,
          typename Encoder = bliss::index::Illumina18QualityScoreCodec<float> >
class QualityScoreGenerationIterator
: public iterator::sliding_window_iterator<BaseIterator,
                                           QualityScoreSlidingWindow<BaseIterator, KMER_SIZE, Encoder > >
{
protected:
  /// The type of the base class
  typedef iterator::sliding_window_iterator<BaseIterator, QualityScoreSlidingWindow<BaseIterator, KMER_SIZE, Encoder > >
    base_class_t;

  /// The difference_type of character offsets
  typedef typename std::iterator_traits<base_class_t>::difference_type diff_type;

  /// The type of the sliding window.
  typedef QualityScoreSlidingWindow<BaseIterator, KMER_SIZE, Encoder > functor_t;

public:
  /// Default constructor.
  QualityScoreGenerationIterator() : base_class_t() {}

  /**
   * @brief   Constructor for the kmer generation iterator using the underlying
   *          iterator `baseBegin` as starting point.
   *
   * @param   baseBegin   An iterator pointing to the first character of the
   *                      sequence to be used for generating k-mers.
   */
  QualityScoreGenerationIterator(const BaseIterator& baseBegin)
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
  QualityScoreGenerationIterator(const BaseIterator& baseBegin, bool initialize_window)
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
  QualityScoreGenerationIterator(const BaseIterator& baseBegin, const functor_t& window)
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
  QualityScoreGenerationIterator(const BaseIterator& baseBegin, const functor_t& window, bool initialize_window)
    : base_class_t(baseBegin, window, initialize_window)  {}
};




} // namespace index
} // namespace bliss

#endif // BLISS_INDEX_QUALITY_SCORE_ITERATOR_HPP
