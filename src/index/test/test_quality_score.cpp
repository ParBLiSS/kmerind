/**
 * @file    test_quality_score.cpp
 * @ingroup
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief
 * @details
 *
 * Copyright (c) 2015 Georgia Institute of Technology.  All Rights Reserved.
 *
 * TODO add License
 */


// include google test
#include <gtest/gtest.h>
//#include <boost/concept_check.hpp>

// include classes to test
#include <index/quality_scores.hpp>
#include <vector>
#include <algorithm>
#include <cassert>

//// Usable AlmostEqual function
//// from http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
//// set maxUlps to be 1M or less.
//bool AlmostEqual2sComplement(float A, float B, int maxUlps)
//{
//    // Make sure maxUlps is non-negative and small enough that the
//    // default NAN won't compare as equal to anything.
//
//    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
//
//    int aInt = *(int*)&A;
//
//    // Make aInt lexicographically ordered as a twos-complement int
//    if (aInt < 0)
//        aInt = 0x80000000 - aInt;
//
//    // Make bInt lexicographically ordered as a twos-complement int
//    int bInt = *(int*)&B;
//
//    if (bInt < 0)
//        bInt = 0x80000000 - bInt;
//
//    int intDiff = abs(aInt - bInt);
//
//    if (intDiff <= maxUlps)
//        return true;
//
//    return false;
//}
//
//bool AlmostEqual2sComplement(double A, double B, int maxUlps)
//{
//    // Make sure maxUlps is non-negative and small enough that the
//    // default NAN won't compare as equal to anything.
//
//    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
//
//    long int aInt = *(long int*)&A;
//
//    // Make aInt lexicographically ordered as a twos-complement int
//    if (aInt < 0)
//        aInt = 0x8000000000000000 - aInt;
//
//    // Make bInt lexicographically ordered as a twos-complement int
//    long int bInt = *(long int*)&B;
//
//    if (bInt < 0)
//        bInt = 0x8000000000000000 - bInt;
//
//    long int intDiff = abs(aInt - bInt);
//
//    if (intDiff <= maxUlps)
//        return true;
//
//    return false;
//}


template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type compare_vectors(const std::vector<T> & first, const std::vector<T> & second) {
  if (first == second) {
    return true;   // if same object
  }
  if (first.size() != second.size()) {
    std::cout << "size not the same" << std::endl;
    return false;  // if size differ
  }

  if (std::equal(first.cbegin(), first.cend(), second.cbegin(), [](T x, T y) {
    return fabs(x - y) < 1.0e-15; }
  )) return true;
  else {
    std::cout << "Vectors not the same" << std::endl;
    for (int i = 0; i < first.size(); ++i) {
      std::cout << (first[i] - second[i]) << std::endl;
    }

    return false;
  }
}

template <typename T>
typename std::enable_if<!std::is_floating_point<T>::value, bool>::type compare_vectors(const std::vector<T> & first, const std::vector<T> & second) {
  if (first == second) return true;   // if same object
  if (first.size() != second.size()) return false;  // if size differ

  if (std::equal(first.cbegin(), first.cend(), second.cbegin())) return true;
  else {

    std::cout << "Vectors not the same" << std::endl;
    for (int i = 0; i < first.size(); ++i) {
      std::cout << first[i] - second[i] << std::endl;
    }

    return false;
  }
}

// templated test function
template<typename CODEC>
void codec_decode(const std::vector<unsigned char> & data, // quality score value
                           std::vector<typename CODEC::value_type>& output) {

  CODEC decoder;

//  std::cout << "decoder LUT SIZE: " << (int)CODEC::size << std::endl;
//  for (int i = 0; i < CODEC::size; ++i) {
//    std::cout << CODEC::lut[i] << ", ";
//  }
//  std::cout << std::endl;

  output.clear();

  for (auto v : data) {
    output.push_back(decoder.decode(v));
  }
}


// templated test function
template<typename CODEC>
void codec_encode(std::vector<typename CODEC::value_type> data, // quality score value
                           std::vector<unsigned char>& output) {

  CODEC decoder;


//  std::cout << "encoder LUT SIZE: " << (int)CODEC::size << std::endl;
//  for (int i = 0; i < CODEC::size; ++i) {
//    std::cout << CODEC::lut[i] << ", ";
//  }
//  std::cout << std::endl;

  output.clear();

  for (auto v : data) {
    output.push_back( decoder.encode(v));
  }
}

template<typename OT, unsigned char MinInput, unsigned char MaxInput, char MinScore>
void direct_decode(const std::vector<unsigned char> & data,
                   std::vector<OT>& output) {

  output.clear();

  OT result;
  for (auto v : data) {
    auto tmp = v - MinInput + MinScore;

    if (tmp == 0) result = std::numeric_limits<OT>::lowest();
    else
      result = std::log2(1.0 - std::exp2(static_cast<OT>(tmp) * log2(10.0) / -10.0));
    output.push_back(result);
  }
}


template<typename OT, unsigned char MinInput, unsigned char MaxInput, char MinScore>
void direct_encode(const std::vector<OT> & data,
                   std::vector<unsigned char>& output) {

  output.clear();

  OT max_log_prob = std::log2(1.0L - std::exp2((static_cast<double>(MaxInput - MinInput + MinScore) - 0.5L) * log2(10.0L) / -10.0L));
  OT min_log_prob = std::log2(1.0L - std::exp2(0.5L * log2(10.0L) / -10.0L));


  unsigned char result;
  for (auto v : data) {
    if (std::isnan(v)) result = MinInput;
    else if (v > max_log_prob) result = MaxInput;
    else if ( v < min_log_prob ) result = MinInput;
    else {
      double tmp = log2(1.0L - exp2(static_cast<double>(v))) * -10.0L / log2(10.0L);
      result = std::round(tmp) + (MinInput - MinScore);
    }
    output.push_back(result);
  }
}


template<typename OT, unsigned char MinInput, unsigned char MaxInput, char MinScore>
void testCodec(const std::string & strdata) {

  std::vector<unsigned char> gold(strdata.begin(), strdata.end());

  // test process.
  std::vector<OT> goldDecoded;
  direct_decode<OT, MinInput, MaxInput, MinScore>(gold, goldDecoded);

  std::vector<unsigned char> goldEncoded;
  direct_encode<OT, MinInput, MaxInput, MinScore>(goldDecoded, goldEncoded);
  bool same = compare_vectors<unsigned char>(goldEncoded, gold);

  if (!same) {
    std::cout << "direct encode/decode: result not same" << std::endl;

    std::cout << std::endl<< "GOLD input: " << std::endl;
    std::copy(gold.begin() , gold.end(), std::ostream_iterator<unsigned char>(std::cout, " "));
    std::cout << std::endl<< "GOLD re-encoded: " << std::endl;
    std::copy(goldDecoded.begin() , goldDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
    std::cout << std::endl<< "GOLD re-encoded: " << std::endl;
    std::copy(goldEncoded.begin() , goldEncoded.end(), std::ostream_iterator<unsigned char>(std::cout, " "));

  }
  EXPECT_TRUE(same);


  // test codec decoding.
  std::vector<OT> codecDecoded;
  codec_decode< bliss::index::QualityScoreCodec<OT, MinInput, MaxInput, MinScore> >(gold, codecDecoded);
  same = compare_vectors<OT>(codecDecoded, goldDecoded);

  if (!same) {
    std::cout << "codec decode: result not same" << std::endl;

    std::cout << std::endl<< "GOLD decoded: " << std::endl;
    std::copy(goldDecoded.begin() , goldDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
    std::cout << std::endl<< "codec decoded: " << std::endl;
    std::copy(codecDecoded.begin() , codecDecoded.end(), std::ostream_iterator<OT>(std::cout, ","));
  }

  EXPECT_TRUE(same);

  std::vector<unsigned char> codecEncodedFromGold;
  codec_encode<bliss::index::QualityScoreCodec<OT, MinInput, MaxInput, MinScore> >(goldDecoded, codecEncodedFromGold);
  same = compare_vectors<unsigned char>(codecEncodedFromGold, gold);

  if (!same) {
    std::cout << "codec encode from goldDecoded : result not same" << std::endl;

    std::cout << std::endl<< "GOLD input: " << std::endl;
    std::copy(gold.begin() , gold.end(), std::ostream_iterator<unsigned char>(std::cout, " "));
    std::cout << std::endl<< "codec encode from goldDecoded: " << std::endl;
    std::copy(codecEncodedFromGold.begin() , codecEncodedFromGold.end(), std::ostream_iterator<unsigned char>(std::cout, " "));

  }
  EXPECT_TRUE(same);

  std::vector<unsigned char> codecEncodedFromCodecDecoded;
  codec_encode<bliss::index::QualityScoreCodec<OT, MinInput, MaxInput, MinScore> >(codecDecoded, codecEncodedFromCodecDecoded);
  same = compare_vectors<unsigned char>(codecEncodedFromCodecDecoded, gold);

  if (!same) {
    std::cout << "codec encode from codecDecoded : result not same" << std::endl;

    std::cout << std::endl<< "GOLD input: " << std::endl;
    std::copy(gold.begin() , gold.end(), std::ostream_iterator<unsigned char>(std::cout, ""));
    std::cout << std::endl<< "codec encode from goldDecoded: " << std::endl;
    std::copy(codecEncodedFromCodecDecoded.begin() , codecEncodedFromCodecDecoded.end(), std::ostream_iterator<unsigned char>(std::cout, " "));

  }
  EXPECT_TRUE(same);

  std::vector<OT> pathologicalValues {
    std::numeric_limits<OT>::quiet_NaN(),
    std::numeric_limits<OT>::signaling_NaN(),
    -std::numeric_limits<OT>::infinity(),
    std::numeric_limits<OT>::lowest(),
    -std::numeric_limits<OT>::denorm_min(),
    -std::numeric_limits<OT>::min(),
    0,
    std::numeric_limits<OT>::denorm_min(),
    std::numeric_limits<OT>::min(),
    std::numeric_limits<OT>::max(),
    std::numeric_limits<OT>::infinity()
  };
  std::vector<unsigned char> pathResults;
//  direct_encode<OT, MinInput, MaxInput, MinScore>(pathologicalValues, pathResults);
  codec_encode<bliss::index::QualityScoreCodec<OT, MinInput, MaxInput, MinScore> >(pathologicalValues, pathResults);
  std::vector<unsigned char> goldPathResults {
    MinInput,
    MinInput,
    MinInput,
    MinInput,
    MaxInput,
    MaxInput,
    MaxInput,
    MaxInput,
    MaxInput,
    MaxInput,
    MaxInput
  };
  same = compare_vectors<unsigned char>(pathResults, goldPathResults);

  if (!same) {
    std::cout << "encode with pathological inputs" << std::endl;

    std::cout << std::endl<< "pathological gold input: " << std::endl;
    std::copy(goldPathResults.begin() , goldPathResults.end(), std::ostream_iterator<OT>(std::cout, ","));
    std::cout << std::endl<< "codec encode pathological input: " << std::endl;
    std::copy(pathResults.begin() , pathResults.end(), std::ostream_iterator<OT>(std::cout, ","));

  }
  EXPECT_TRUE(same);

}


/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGeneration, TestIllunima18QualityScoreDouble)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testCodec<double, 33, 126, 0>(strdata);
}




/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGeneration, TestIllunima13QualityScoreDouble)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testCodec<double, 64, 126, 0>(strdata);
}



/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGeneration, TestIllunima15QualityScoreDouble)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "CDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testCodec<double, 67, 126, 3>(strdata);

}


/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGeneration, TestIllunima18QualityScoreFloat)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testCodec<float, 33, 126, 0>(strdata);

}




/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGeneration, TestIllunima13QualityScoreFloat)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testCodec<float, 64, 126, 0>(strdata);

}



/**
 * Test k-mer generation with 2 bits for each character
 */
TEST(QualityScoreGeneration, TestIllunima15QualityScoreFloat)
{
  // test sequence: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  std::string strdata = "CDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  testCodec<float, 67, 126, 3>(strdata);

}

