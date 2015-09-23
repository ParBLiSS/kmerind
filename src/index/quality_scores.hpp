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
#include <algorithm>

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
 *  supports float or double as output format.  float representation of probability starts to fail at Q of 75 in [0, 96)
 *  in Sangar encoding of raw base reads (0 to 93).  For kmer composite scores, this values needs to be more precise
 *  so double is appropriate.
 *
 * Creates a look up table instead of computing (by calling std::log, exp, etc.).  Most of usages is to look up
 * phred ascii score to log probability of base being correct.
 * Using LUT makes faster by about 3 seconds on a 35MB read file.
 *
 *  to avoid the use of constexpr versions of log2, round, and exp2, we can pre-compute the entire phred score table.
 * This particular version is for PhredScore only, since Solexa/Illumina1.0 uses odds instead of probability for quality score
 * and the values are NOT compatible.
 *
 *
 * The input scores range between 0 and 93, but most input values will be between 0 and 41.
 * For composite scores (of a kmer, for example), we could have higher scores.  To bound the LUT size, we choose
 * 96 entries, which is 12 cacheline's of doubles.  mapping to ASCII is restricted to 0 to 93, corresponding to ASCII 33 to 126
 *
 *
 *  2 tables are created.  first is for decoding (from ASCII), and second is for encoding (to ASCII).
 *  In both there is assumption that a range of probabilities (i.e. log2(p_corr)) are mapped to the same character via ROUNDING (instead of CEIL or FLOOR)
 *  decoding table maps phred score to log2(p_correct) at integer q values.
 *  encoding table contains phred scores at boundary of rounding, i.e. val - std::numeric_limits<double>::round_error, for lookup.
 *
 *  table will store values log2(prob_correct). this is mathematically log2(1.0 - p_err) = log2(1.0 - 10^(q/(-10)))
 *    = log2(1.0 - exp2(log2(10) * (q / (-10)))) = log2(1.0 - exp2(q * (log2(10) / (-10))));
 *  NOTE: q = 0 results in log2(p_correct) = -inf.  change to log2(p_correct) = std::numeric_limits<double>::lowest();
 *
 *  operator() is no longer needed as the lookup table is precomputed and hardcoded.
 *  encode() and decode() methods are used at runtime.
 *
 *  Decode(q_score) logic computes the appropriate LUT index based on MinInput, MaxInput, and MinScore, then look-up.
 *  Encode(log2(p_correct)) logic performs a binary search in lookup table to find the first >= (same as lowerbound) entry, after which
 *    the index is transformed based on MinInput, MaxInput, and MinScore back to an ASCII character.
 *
 */
template<typename OutT, unsigned char MinInput, unsigned char MaxInput, char MinScore >
struct QualityScoreCodec
{
private:
    // OutT has to be float or double
    static_assert(std::is_floating_point<OutT>::value, "Quality Score Codec output needs to be floating point type");
    static_assert(MaxInput >= MinInput, "Quality Score Input range is invalid");
    static_assert(MinScore >= 0, "Minimum score should at least 0 for Phred Scores");

public:
    /// number of possible quality-score character values
    static constexpr unsigned char size = MaxInput - MinInput + 1;

protected:
    /// Type of the lookup-table
    typedef std::array<OutT, 96> LUTType;

public:
    static constexpr OutT min_log_prob = std::numeric_limits<OutT>::lowest();

    /*
     * @brief a lookup table (static constexpr array) for ASCII to log2(probability(base_correct))
     * @details Performance of lookup-table: saves about 3 seconds on a 35MB read file.
     */
    // generated from  std::log2(1.0L - std::exp2(static_cast<long double>(q) * std::log2(10.0L) / (-10.0L)))
    // (use it in a loop in a c++ program)
    static constexpr LUTType DecodeLUT = {{
       // log2(p_correct)                     // q
       std::numeric_limits<OutT>::lowest(), // 0  value is instead of -inf
       MinScore > 1 ? std::numeric_limits<OutT>::lowest() : -2.28158434133843178,                  // 1
       MinScore > 2 ? std::numeric_limits<OutT>::lowest() : -1.43814051613477932,                  // 2
       -1.00342970560804731,                  // 3
       -0.73242146536126580,                  // 4
       -0.54841225460816377,                  // 5
       -0.41732577916269975,                  // 6
       -0.32107396843632689,                  // 7
       -0.24894651201558545,                  // 8
       -0.19411744585050233,                  // 9
       -0.15200309344504998,                  // 10
       -0.11940509171847160,                  // 11
       -0.09402645646885732,                  // 12
       -0.07418088913876119,                  // 13
       -0.05860926129544923,                  // 14
       -0.04635894788914435,                  // 15
       -0.03670176875519015,                  // 16
       -0.02907660209329712,                  // 17
       -0.02304830733610452,                  // 18
       -0.01827774903217653,                  // 19
       -0.01449956969511508,                  // 20
       -0.01150549046750825,                  // 21
       -0.00913162905170873,                  // 22
       -0.00724878356560029,                  // 23
       -0.00575493542826523,                  // 24
       -0.00456943101694318,                  // 25
       -0.00362844512928738,                  // 26
       -0.00288143060891852,                  // 27
       -0.00228833140892072,                  // 28
       -0.00181738966764433,                  // 29
       -0.00144341686966872,                  // 30
       -0.00114642878575454,                  // 31
       -0.00091056632636779,                  // 32
       -0.00072324159138011,                  // 33
       -0.00057446159692902,                  // 34
       -0.00045629237978688,                  // 35
       -0.00036243413137627,                  // 36
       -0.00028788422589034,                  // 37
       -0.00022866987625544,                  // 38
       -0.00018163597839181,                  // 39
       -0.00014427671804504,                  // 40
       -0.00011460189214367,                  // 41
       -0.00009103077504641,                  // 42
       -0.00007230784565469,                  // 43
       -0.00005743586735910,                  // 44
       -0.00004562274434512,                  // 45
       -0.00003623931612729,                  // 46
       -0.00002878583764822,                  // 47
       -0.00002286535668724,                  // 48
       -0.00001816256881128,                  // 49
       -0.00001442702254412,                  // 50
       -0.00001145977956494,                  // 51
       -0.00000910281903642,                  // 52
       -0.00000723062148460,                  // 53
       -0.00000574348383962,                  // 54
       -0.00000456220951173,                  // 55
       -0.00000362389064942,                  // 56
       -0.00000287855791882,                  // 57
       -0.00000228651936105,                  // 58
       -0.00000181624659170,                  // 59
       -0.00000144269576224,                  // 60
       -0.00000114597386021,                  // 61
       -0.00000091027931907,                  // 62
       -0.00000072306051771,                  // 63
       -0.00000057434735503,                  // 64
       -0.00000045622030196,                  // 65
       -0.00000036238865532,                  // 66
       -0.00000028785553343,                  // 67
       -0.00000022865177303,                  // 68
       -0.00000018162455628,                  // 69
       -0.00000014426951130,                  // 70
       -0.00000011459734506,                  // 71
       -0.00000009102790606,                  // 72
       -0.00000007230603546,                  // 73
       -0.00000005743472521,                  // 74
       -0.00000004562202370,                  // 75
       -0.00000003623886144,                  // 76
       -0.00000002878555076,                  // 77
       -0.00000002286517567,                  // 78
       -0.00000001816245460,                  // 79
       -0.00000001442695048,                  // 80
       -0.00000001145973410,                  // 81
       -0.00000000910279035,                  // 82
       -0.00000000723060338,                  // 83
       -0.00000000574347242,                  // 84
       -0.00000000456220231,                  // 85
       -0.00000000362388610,                  // 86
       -0.00000000287855505,                  // 87
       -0.00000000228651755,                  // 88
       -0.00000000181624545,                  // 89
       -0.00000000144269504,                  // 90
       -0.00000000114597341,                  // 91
       -0.00000000091027903,                  // 92
       -0.00000000072306034,                  // 93
       0.0,                                   // 94.  no error. log2(1) = 0.  these values are used in accumulation, so 0.0 is appropriate.
       0.0                                    // 95.  no error. log2(1) = 0
    }};

    // generated from std::log2(1.0L - std::exp2((static_cast<long double>(q) - std::numeric_limits<long double>::round_error()) * std::log2(10.0L) / (-10.0L)))
    // (use it in a loop in a c++ program)
    // note that this does not address rounding errors (e.g. some values correspond to X.4999999999).  hopefully the errors occurs infrequently.
    static constexpr LUTType EncodeLUT = {{
       // log2(p_correct)                     // q
       std::numeric_limits<OutT>::lowest(),   // 0   '!'   // instead of -inf
       MinScore > 1 ? std::numeric_limits<OutT>::lowest() : -3.20092513941326206,                  // 1   '"'
       MinScore > 2 ? std::numeric_limits<OutT>::lowest() : -1.77569188557725199,                  // 2   '#'
       -1.19212192855224533,                  // 3   '$'
       -0.85382338930387750,                  // 4   '%'
       -0.63221159560848881,                  // 5   '&'
       -0.47761936582863019,                  // 6   '''
       -0.36563370342444848,                  // 7   '('
       -0.28248775124226952,                  // 8   ')'
       -0.21969620844535423,                  // 9   '*'
       -0.17169638509320007,                  // 10  '+'
       -0.13467515789197097,                  // 11  ','
       -0.10593052313711612,                  // 12  '-'
       -0.08349909941006277,                  // 13  '.'
       -0.06592644712255149,                  // 14  '/'
       -0.05211894372310457,                  // 15  '0'
       -0.04124465628238184,                  // 16  '1'
       -0.03266493703250296,                  // 17  '2'
       -0.02588599985503053,                  // 18  '3'
       -0.02052390673788567,                  // 19  '4'
       -0.01627880192555448,                  // 20  '5'
       -0.01291567449858786,                  // 21  '6'
       -0.01024982328062335,                  // 22  '7'
       -0.00813576730234849,                  // 23  '8'
       -0.00645871781187414,                  // 24  '9'
       -0.00512797793639086,                  // 25  ':'
       -0.00407180773406806,                  // 26  ';'
       -0.00323341259107256,                  // 27  '<'
       -0.00256779869522924,                  // 28  '='
       -0.00203930154501344,                  // 29  '>'
       -0.00161963926239068,                  // 30  '?'
       -0.00128637663601342,                  // 31  '@'
       -0.00102171157323732,                  // 32  'A'
       -0.00081151523680720,                  // 33  'B'
       -0.00064457217386949,                  // 34  'C'
       -0.00051197835084737,                  // 35  'D'
       -0.00040666401726464,                  // 36  'E'
       -0.00032301534686808,                  // 37  'F'
       -0.00025657430233395,                  // 38  'G'
       -0.00020380048520041,                  // 39  'H'
       -0.00016188212788732,                  // 40  'I'
       -0.00012858606107082,                  // 41  'J'
       -0.00010213860270725,                  // 42  'K'
       -0.00008113098530057,                  // 43  'L'
       -0.00006444425964465,                  // 44  'M'
       -0.00005118965985303,                  // 45  'N'
       -0.00004066124378008,                  // 46  'O'
       -0.00003229828038102,                  // 47  'P'
       -0.00002565537697460,                  // 48  'Q'
       -0.00002037875303607,                  // 49  'R'
       -0.00001618739541103,                  // 50  'S'
       -0.00001285809038534,                  // 51  'T'
       -0.00001021353487668,                  // 52  'U'
       -0.00000811289322249,                  // 53  'V'
       -0.00000644429642524,                  // 54  'W'
       -0.00000511888425213,                  // 55  'X'
       -0.00000406607280814,                  // 56  'Y'
       -0.00000322979549985,                  // 57  'Z'
       -0.00000256551716728,                  // 58  '['
       -0.00000203786234998,                  // 59  '\'
       -0.00000161873136793,                  // 60  ']'
       -0.00000128580388162,                  // 61  '^'
       -0.00000102135023388,                  // 62  '_'
       -0.00000081128726925,                  // 63  '`'
       -0.00000064442834717,                  // 64  'a'
       -0.00000051188760790,                  // 65  'b'
       -0.00000040660676512,                  // 66  'c'
       -0.00000032297922461,                  // 67  'd'
       -0.00000025655151143,                  // 68  'e'
       -0.00000020378610546,                  // 69  'f'
       -0.00000016187305506,                  // 70  'g'
       -0.00000012858033659,                  // 71  'h'
       -0.00000010213499085,                  // 72  'i'
       -0.00000008112870640,                  // 73  'j'
       -0.00000006444282176,                  // 74  'k'
       -0.00000005118875262,                  // 75  'l'
       -0.00000004066067136,                  // 76  'm'
       -0.00000003229791921,                  // 77  'n'
       -0.00000002565514909,                  // 78  'o'
       -0.00000002037860925,                  // 79  'p'
       -0.00000001618730469,                  // 80  'q'
       -0.00000001285803314,                  // 81  'r'
       -0.00000001021349876,                  // 82  's'
       -0.00000000811287043,                  // 83  't'
       -0.00000000644428205,                  // 84  'u'
       -0.00000000511887518,                  // 85  'v'
       -0.00000000406606708,                  // 86  'w'
       -0.00000000322979189,                  // 87  'x'
       -0.00000000256551489,                  // 88  'y'
       -0.00000000203786091,                  // 89  'z'
       -0.00000000161873046,                  // 90  '{'
       -0.00000000128580331,                  // 91  '|'
       -0.00000000102134987,                  // 92  '}'
       -0.00000000081128704,                  // 93  '~'
       std::numeric_limits<OutT>::max(),    // 94       want to force to choose 93 during search
       std::numeric_limits<OutT>::max()     // 95       want to force to choose 93 during search
    }};



    /// Type of the returned log-probability score (float or double)
    typedef OutT value_type;


    /**
     * @brief Returns the log_2 probability that the base is correct converted
     *        from the given Phred Score.
     *
     * @param score The raw, encoded score as ASCII character.
     *
     * @return The log2 probability that the base with the given quality score
     *         is correct.
     */
    inline static OutT decode(const unsigned char score)
    {
      return DecodeLUT[score - MinInput];  // less than MinScore  - DecodeLUT values are same as idx = 0;
    }

    /**
     * @brief Returns the Sanger/Phred score for the given log2 probability
     *        of a base/k-mer being correct.
     *
     * @param log2_prob_correct The log2 prob(correct).
     * @return The Sanger/Phred score corresponding to the given probability.
     */
    inline static unsigned char encode(const OutT log2_prob_correct)
    {
    	if (std::isnan(log2_prob_correct)) return MinScore == 0 ? MinInput : (MinInput + MinScore - 1);
    	if (std::isinf(log2_prob_correct)) {
    		if (std::signbit(log2_prob_correct))
    			return MinScore == 0 ? MinInput : (MinInput + MinScore - 1);
    		else return MaxInput;
    	}
    	if (log2_prob_correct == std::numeric_limits<OutT>::lowest())
    		return MinScore == 0 ? MinInput : (MinInput + MinScore - 1);
    	if (log2_prob_correct == std::numeric_limits<OutT>::max())
    		return MaxInput;


      // do binary search (search for upper bound:  *it > val)
      auto it = std::upper_bound(EncodeLUT.cbegin(), EncodeLUT.cend(), log2_prob_correct);

      if (it == EncodeLUT.cbegin()) {
        // at beginning, no entry in LUT is <= val.  assign lowest quality score char
        return MinScore == 0 ? MinInput : (MinInput + MinScore - 1);

      } else if (it == EncodeLUT.cend()) {
        // everything's less, so set to max.
        return MaxInput;

      } else {
        // found it.  go back 1 since it is at the upper bound (exclusive). (since not at cbegin(), it's okay)
        return std::min<unsigned char>(MaxInput, MinInput + std::distance(EncodeLUT.cbegin(), it) - 1);
      }
    }
};


///**
// * @class QualityScoreCodec, converts between Phred Scores (encoded in ASCII) to Log(Pr(correct bases)).
// * @details  see http://en.wikipedia.org/wiki/FASTQ_format.
// *
// *  supports float or double as output format.  float representation of probability starts to fail at Q of 98
// *  in Sangar encoding of raw base reads (0 to 93).  For kmer composite scores, this values needs to be more precise
// *
// * Creates a look up table instead of computing (by calling std::log, exp, etc.).  Most of usages is to look up
// * phred ascii score to log probability of base being correct.
// * Using LUT makes faster by about 3 seconds on a 35MB read file.
// *
// *  to avoid the use of constexpr versions of log2, round, and exp2, we can pre-compute the entire encode and decode tables.
// *
// *
// *
// *  we store log2(prob_correct), and index has range [MinScore, MaxInput - MinInput + MinScore]
// *           valid phred score range is [0, MaxInput - MinInput + MinScore]
// *           prob_correct has range [0, 1], so the stored values then have range [-inf, 0].
// *           given phred score is discrete, we can say that phred score range that requires actual computation
// *              is [0.5, MaxInput - MinInput + MinScore - 0.5]; outside of this range we can do lookup, i.e
// *              phred score < 0.5 => MinLog2, and phred_score > MaxInput-MinInput +MinScore - 0.5 => MaxLog2
// *
// *          computing log2(probcorrect) from phred score is straight forward.
// *          computing phred score from log2(probcorrect) can benefit from tighter bounds to avoid numeric errors
// *          and boundary cases.
// *
// *          log2(probcorrect) on the side of -inf of range maps to maximum phred score.  given max phred score,
// *          we can then compute the maximum log2(prob correct)
// *             max phred score = MaxInput - MinInput + MinScore - 0.5
// *             so max log2(prob correct) =
// *          at the -inf end, math is okay.  at 0 end, we can check to see if input is > x, then MaxInput
// *          x is derived from phred score of 0.5, where p = 10^(-0.05), so corresponding input is
// *          log_2(1-10^(-0.05)) = -3.200925139, so check for input < -3.200915139, to set to MinInput
// */
//template<typename OutT, unsigned char MinInput, unsigned char MaxInput, char MinScore >
//struct QualityScoreCodec
//{
//private:
//    // OutT has to be float or double
//    static_assert(std::is_floating_point<OutT>::value, "Quality Score Codec output needs to be floating point type");
//    static_assert(MaxInput >= MinInput, "Quality Score Input range is invalid");
//
//
//public:
//    /// number of possible quality-score character values
//    static constexpr unsigned char size = MaxInput - MinInput + 1;
//
//protected:
//    /// Type of the lookup-table
//    typedef std::array<OutT, size> LUTType;
//
//    /// Conversion coefficient for Phred-Score -> probablity(base_incorrect)
//    static constexpr double log2_10DivNeg10 = -0.33219280948873623;
//    //static constexpr double log2_10DivNeg10 = std::log2(10.0L) / -10.0L;
//    /// Conversion coefficient for probablity(base_incorrect) -> Phred-Score
//    static constexpr double neg10DivLog2_10 = -3.01029995663981195;
//    //static constexpr double neg10DivLog2_10 = -10.0L / std::log2(10.0L);
//
//
//public:
//    /*
//     * @brief a lookup table (static constexpr array) for ASCII to log2(probability(base_correct))
//     * @details Performance of lookup-table: saves about 3 seconds on a 35MB read file.
//     */
//    static constexpr LUTType lut = bliss::utils::make_array<size>(QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>());
//
//
//    /// Type of the returned log-probability score (float or double)
//    typedef OutT value_type;
//
//    /// maximum possible log probability value given parameters
//    static constexpr OutT max_log_prob = std::log2(1.0L - std::exp2((static_cast<double>(MaxInput - MinInput + MinScore) - 0.5L) * log2_10DivNeg10));
//    /// minimum possible log probability value given parameters
//    static constexpr OutT min_log_prob = -3.20092513941326206;
//    // static constexpr OutT min_log_prob = std::log2(1.0L - std::exp2(0.5L * log2_10DivNeg10));
//    /// log probability value for a quality score that show error in base read.
//    static constexpr OutT error_log_prob = std::numeric_limits<OutT>::lowest();
//
//    /**
//     * @brief   Return the log2 probability of correctness for the provided Phred score value.
//     * @details  used by bliss::utils::make_array to create LUT.
//     *      requires input be size_t type, as this transforms from 0..N to log prob of correct base.
//     *      if the quality score is not a phred score, then template specialization is necesssary.
//     *
//     *        idx is offset from first array entry; idx = phred - MinScore
//     */
//    constexpr OutT operator()(const size_t &idx) const
//    {
//      // needs to be public for bliss::utils::make_array
//
//      // some limits: v / -10 has to be negative as this becomes probability, so idx > 0
//      return
//          // if the score is invalid, assume the base is incorrect, set the value to  -inf
//          static_cast<double>(idx) <= - static_cast<double>(MinScore) ? error_log_prob :  // min score less than or equal to 0 results in invalid value.
//          idx >= size ? 0.0 :  // maximum phred score means prob of correct is 1., log of which is 0.
//          // Phred Score to log2(prob(correct)) per base:
//              // log_2(p) / log_2(10) = log_10(p)
//          std::log2(1.0L - std::exp2((static_cast<double>(idx) + static_cast<double>(MinScore)) * log2_10DivNeg10));
//    }
//
//    /**
//     * @brief Returns the log_2 probability that the base is correct converted
//     *        from the given Phred Score.
//     *
//     * @param score The raw, encoded score as ASCII character.
//     *
//     * @return The log2 probability that the base with the given quality score
//     *         is correct.
//     */
//    inline static constexpr OutT decode(const unsigned char score)
//    {
//        return lut[score - MinInput];  // for Illumina 1.0 score, MinScore would be -5.
//    }
//
//    /**
//     * @brief Returns the Sanger/Phred score for the given log2 probability
//     *        of a base/k-mer being correct.
//     *
//     * @param log2_prob_correct The log2 prob(correct).
//     * @return The Sanger/Phred score corresponding to the given probability.
//     */
//    inline static constexpr unsigned char encode(const OutT log2_prob_correct)
//    {
//        return
//            std::isnan(log2_prob_correct) ? MinInput :
//            log2_prob_correct < min_log_prob ? MinInput :  // correspond to phred of 0.5
//            log2_prob_correct > max_log_prob ? MaxInput :  // correspond to phred of size - 0.5
//            MinInput +
//                std::min<char>(size - 1,
//                               std::max<char>(0,
//                                              static_cast<char>(std::round(neg10DivLog2_10 *
//                                                                           std::log2(1.0L - std::exp2(static_cast<double>(log2_prob_correct))) -
//                                                                           static_cast<OutT>(MinScore)))
//                                              )
//                               );
//    }
//};




// NOTE:  Define quality score lookup table here.  DO NOT define in separate cpp file or where it mightbe used.
/// define (not declare or initialize) the quality score lookup table.
template<typename OutT, unsigned char MinInput, unsigned char MaxInput, char MinScore>
constexpr typename QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>::LUTType QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>::DecodeLUT;

/// define (not declare or initialize) the quality score lookup table.
template<typename OutT, unsigned char MinInput, unsigned char MaxInput, char MinScore>
constexpr typename QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>::LUTType QualityScoreCodec<OutT, MinInput, MaxInput, MinScore>::EncodeLUT;


/// Illumina 1.8 quality score converstion presets.  convenience typedef.
template<typename OutT>
using Illumina18QualityScoreCodec = QualityScoreCodec<OutT, 33, 126, 0>;

/// Sanger quality score converstion presets. same as Illumina 1.8 convenience typedef.
template<typename OutT>
using SangerQualityScoreCodec = QualityScoreCodec<OutT, 33, 126, 0>;


/// Illumina 1.3 quality score converstion presets.  convenience typedef
template<typename OutT>
using Illumina13QualityScoreCodec = QualityScoreCodec<OutT, 64, 126, 0>;

/// Illumina 1.5 quality score converstion presets.  convenience typedef
template<typename OutT>
using Illumina15QualityScoreCodec = QualityScoreCodec<OutT, 64, 126, 3>;  // special Q val of 2 indicate bases should not be used, not a quality measurement

// illumina 1.0/Solexa uses a different quality score measurement (-10 log_10(p/1-p)), not phred score.  for now, disable.
//template<typename OutT>
//using SolexaQualityScoreCodec = QualityScoreCodec<OutT, 59, 126, -5>;



} // namespace index
} // namespace bliss

#endif // BLISS_INDEX_QUALITY_SCORES_HPP
