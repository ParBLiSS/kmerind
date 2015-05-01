/**
 * @file    quality_scores.hpp
 * @ingroup index
 * @author  Tony Pan <tpan7@gatech.edu>
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements common functions for handleing of quality scores.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 *
 * TODO add Licence
 */

#ifndef BLISS_INDEX_QUALITY_SCORES_HPP
#define BLISS_INDEX_QUALITY_SCORES_HPP

#include <cstdlib>

#include <numeric>
#include <cmath>
#include <array>
#include <limits>

#include "utils/constexpr_array.hpp"


namespace bliss
{
namespace index
{
template <typename FileFormat>
struct QualityType;


/**
 * @class QualityScoreCodec, converts between Phred Scores (encoded in ASCII) to Log(Pr(correct bases)).
 * @details  see http://en.wikipedia.org/wiki/FASTQ_format.
 *
 *  supports float or double as output format.  float representation of probability starts to fail at Q of 98
 *  in Sangar encoding of raw base reads (0 to 93).  For kmer composite scores, this values needs to be more precise
 *
 * Creates a look up table instead of computing (by calling std::log, exp, etc.).  Most of usages is to look up
 * phred ascii score to log probability of base being correct.
 * Using LUT makes faster by about 3 seconds on a 35MB read file.
 *
 *  we store log2(prob_correct), and index has range [MinScore, MaxInput - MinInput + MinScore]
 *           valid phred score range is [0, MaxInput - MinInput + MinScore]
 *           prob_correct has range [0, 1], so the stored values then have range [-inf, 0].
 *           given phred score is discrete, we can say that phred score range that requires actual computation
 *              is [0.5, MaxInput - MinInput + MinScore - 0.5]; outside of this range we can do lookup, i.e
 *              phred score < 0.5 => MinLog2, and phred_score > MaxInput-MinInput +MinScore - 0.5 => MaxLog2
 *
 *          computing log2(probcorrect) from phred score is straight forward.
 *          computing phred score from log2(probcorrect) can benefit from tighter bounds to avoid numeric errors
 *          and boundary cases.
 *
 *          log2(probcorrect) on the side of -inf of range maps to maximum phred score.  given max phred score,
 *          we can then compute the maximum log2(prob correct)
 *             max phred score = MaxInput - MinInput + MinScore - 0.5
 *             so max log2(prob correct) =
 *          at the -inf end, math is okay.  at 0 end, we can check to see if input is > x, then MaxInput
 *          x is derived from phred score of 0.5, where p = 10^(-0.05), so corresponding input is
 *          log_2(1-10^(-0.05)) = -3.200925139, so check for input < -3.200915139, to set to MinInput
 */
template<typename OutT, unsigned char MinInput, unsigned char MaxInput, char MinScore >
struct QualityScoreCodec
{
private:
    // OutT has to be float or double
    static_assert(std::is_floating_point<OutT>::value, "Quality Score Codec output needs to be floating point type");
    static_assert(MaxInput >= MinInput, "Quality Score Input range is invalid");


public:
    /// number of possible quality-score character values
    static constexpr unsigned char size = MaxInput - MinInput + 1;

protected:
    /// Type of the lookup-table
    typedef std::array<OutT, size> LUTType;

    /// Conversion coefficient for Phred-Score -> probablity(base_incorrect)
    static constexpr double log2_10DivNeg10 = std::log2(10.0L) / -10.0L;
    /// Conversion coefficient for probablity(base_incorrect) -> Phred-Score
    static constexpr double neg10DivLog2_10 = -10.0L / std::log2(10.0L);


public:
    /*
     * @brief a lookup table (static constexpr array) for ASCII to log2(probability(base_correct))
     * @details Performance of lookup-table: saves about 3 seconds on a 35MB read file.
     */
    static constexpr LUTType lut = bliss::utils::make_array<size>(QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>());


    /// Type of the returned log-probability score (float or double)
    typedef OutT value_type;

    /// maximum possible log probability value given parameters
    static constexpr OutT max_log_prob = std::log2(1.0L - std::exp2((static_cast<double>(MaxInput - MinInput + MinScore) - 0.5L) * log2_10DivNeg10));
    /// minimum possible log probability value given parameters
    static constexpr OutT min_log_prob = std::log2(1.0L - std::exp2(0.5L * log2_10DivNeg10));
    /// log probability value for a quality score that show error in base read.
    static constexpr OutT error_log_prob = std::numeric_limits<OutT>::lowest();

    /**
     * @brief   Return the log2 probability of correctness for the provided Phred score value.
     * @details  used by bliss::utils::make_array to create LUT.
     * 			requires input be size_t type, as this transforms from 0..N to log prob of correct base.
     * 			if the quality score is not a phred score, then template specialization is necesssary.
     *
     * 		    idx is offset from first array entry; idx = phred - MinScore
     */
    constexpr OutT operator()(const size_t &idx) const
    {
      // needs to be public for bliss::utils::make_array

      // some limits: v / -10 has to be negative as this becomes probability, so idx > 0
      return
          // if the score is invalid, assume the base is incorrect, set the value to  -inf
          static_cast<double>(idx) <= - static_cast<double>(MinScore) ? error_log_prob :  // min score less than or equal to 0 results in invalid value.
          idx >= size ? 0.0 :  // maximum phred score means prob of correct is 1., log of which is 0.
          // Phred Score to log2(prob(correct)) per base:
              // log_2(p) / log_2(10) = log_10(p)
          std::log2(1.0L - std::exp2((static_cast<double>(idx) + static_cast<double>(MinScore)) * log2_10DivNeg10));
    }

    /**
     * @brief Returns the log_2 probability that the base is correct converted
     *        from the given Phred Score.
     *
     * @param score The raw, encoded score as ASCII character.
     *
     * @return The log2 probability that the base with the given quality score
     *         is correct.
     */
    inline static constexpr OutT decode(const unsigned char score)
    {
        return lut[score - MinInput];  // for Illumina 1.0 score, MinScore would be -5.
    }

    /**
     * @brief Returns the Sanger/Phred score for the given log2 probability
     *        of a base/k-mer being correct.
     *
     * @param log2_prob_correct The log2 prob(correct).
     * @return The Sanger/Phred score corresponding to the given probability.
     */
    inline static constexpr unsigned char encode(const OutT log2_prob_correct)
    {
        return
            std::isnan(log2_prob_correct) ? MinInput :
            log2_prob_correct < min_log_prob ? MinInput :  // correspond to phred of 0.5
            log2_prob_correct > max_log_prob ? MaxInput :  // correspond to phred of size - 0.5
            MinInput +
                std::min<char>(size - 1,
                std::max<char>(0,
                       static_cast<char>(std::round(neg10DivLog2_10 *
                                                    std::log2(1.0L -
                                                              std::exp2(static_cast<double>(log2_prob_correct))) -
                                                              static_cast<OutT>(MinScore))
                                         )
                                     )
                                     );
    }
};

// NOTE:  Do it here is fine.  do it in class file or even where it's used is
// NOT.  probably due to templating.
/// define (not declare or initialize) the quality score lookup table.
template<typename OutT, unsigned char MinInput, unsigned char MaxInput, char MinScore>
constexpr typename QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>::LUTType QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>::lut;

/// Illumina 1.8 quality score converstion presets.  convenience typedef.
template<typename OutT>
using Illumina18QualityScoreCodec = QualityScoreCodec<OutT, 33, 126, 0>;

/// Illumina 1.3 quality score converstion presets.  convenience typedef
template<typename OutT>
using Illumina13QualityScoreCodec = QualityScoreCodec<OutT, 64, 126, 0>;

/// Illumina 1.5 quality score converstion presets.  convenience typedef
template<typename OutT>
using Illumina15QualityScoreCodec = QualityScoreCodec<OutT, 67, 126, 3>;  // ignore special Q val of 2.

// illumina 1.0/Solexa uses a different quality score measurement, not phred score.  for now, disable.
//template<typename OutT>
//using Illumina10QualityScoreCodec = QualityScoreCodec<OutT, 64, 126, 0>;



} // namespace index
} // namespace bliss

#endif // BLISS_INDEX_QUALITY_SCORES_HPP
